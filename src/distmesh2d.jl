########################################################################
# StaticVector shorthands

const SF2 = SVector{2,Float64}
const SI2 = SVector{2,Int32}
const SI3 = SVector{3,Int32}

########################################################################
# Utility functions

extract_elements(tri::Delaunator.Triangulation{I}) where {I} =
    reinterpret(SVector{3, I}, triangles(tri))
delaunay(p) = extract_elements(triangulate(p))

barvectors(p, bars) = [ p[bar[1]] - p[bar[2]] for bar in bars ]

function init_nodes(bbox, h0)
    # Create initial distribution in bounding box (equilateral triangles)
    xx = bbox[1][1]:h0:bbox[2][1]
    yy = bbox[1][2]:h0*√3/2:bbox[2][2]
    p = [ SF2(x + iseven(iy) * h0/2, y) for x in xx for (iy,y) in enumerate(yy) ]
end

function nodes_rejection(p, pfix, dfcn, hfcn, geps)
    # Remove nodes outside the region, apply the rejection method
    p = p[dfcn.(p) .< geps]
    r0 = [ 1 / hfcn(pp)^2 for pp in p ]
    p = p[rand(length(r0)) .< r0./maximum(r0)]
    p = vcat(pfix, setdiff(p, pfix))
end

function retriangulate(p, dfcn, geps)
    # Retriangulation by the Delaunay algorithm, remove outside triangles
    t = delaunay(p)
    t = filter(tt -> dfcn(sum(p[tt]) / 3) < -geps, t)
    # Find all bars as the unique triangle edges
    bars = [ tt[ix] for tt in t for ix in SI2[(1,2),(2,3),(3,1)] ]
    bars = unique(sort.(bars))
    t, bars
end

function desiredlengths(p, bars, L, hfcn, Fscale)
    # Scaled and normalized desired bar lengths
    barmid = [ (p[bar[1]] + p[bar[2]]) / 2 for bar in bars ]
    hbars = hfcn.(barmid)
    L0 = hbars * Fscale * sqrt(sum(L .^ 2) / sum(hbars .^ 2))
end

function density_control(p, L, L0, bars, nfix)
    # Density control - remove nodes that are too close
    pkeep = trues(length(p))
    foreach(bar -> pkeep[bar] .= false, bars[L0[:].>2*L[:]])
    pkeep[1:nfix] .= true
    p = p[pkeep]
end

function total_node_forces(L, L0, barvec, bars, np, nfix)
    # Find all bar forces and accumulate at all nodes
    F = max.(L0 .- L, 0.0)
    Fvec = F ./ L .* barvec
    Ftot = [ SF2(0.0,0.0) for ip = 1:np ]
    for ibar in eachindex(bars)
        Ftot[bars[ibar][1]] += Fvec[ibar]
        Ftot[bars[ibar][2]] -= Fvec[ibar]
    end
    Ftot[1:nfix] .*= 0.0
    Ftot
end

function project_nodes!(p, dfcn, deps)
    d = dfcn.(p)
    ix = findall(d .> 0)
    numgrad(f,x,fx) = (SF2(f(x + SF2(deps,0)), f(x + SF2(0,deps))) .- fx) / deps
    dgrad = [ numgrad(dfcn, p[i], d[i]) for i in ix ]
    @. p[ix] -= d[ix] * dgrad / norm(dgrad)^2
    d
end

########################################################################
# Main function

"""
        distmesh2d(dfcn, hfcn, h0, bbox, pfix=[]; plotting=false)

2-D Mesh Generator using Distance Functions.

        `p`: Node positions (Vector of (x,y) coordinates)
        't': Triangle indices (Vector of (i1,i2,i3) indices)
     'dfcn': Distance function d(p)
     'hfcn': Scaled edge length function h(p)
       'h0': Initial edge length
     'bbox': Bounding box ((xmin,ymin), (xmax,ymax))
     'pfix': Fixed node positions (Vector of (x,y) coordinates)

# Examples

### Example: (Uniform Mesh on Unit Circle)
```julia
dfcn(p) = sqrt(sum(p.^2)) - 1
p,t = distmesh2d(dfcn, huniform, 0.2, ((-1,-1), (1,1)), plotting=true);
```
#### Example: (Rectangle with circular hole, refined at circle boundary)
```julia
dfcn2(p) = ddiff(drectangle(p, -1, 1, -1, 1), dcircle(p, r=0.5))
hfcn2(p) = 0.05 + 0.3 * dcircle(p, r=0.5);
bbox = ((-1,-1), (1,1))
pfix = [(-1,-1), (-1,1), (1,-1), (1,1)]
p,t = distmesh2d(dfcn2, hfcn2, 0.05, bbox, pfix, plotting=true);
```
"""
function distmesh2d(dfcn, hfcn, h0, bbox, pfix=SF2[];
                    plotting=false,          # Optional live plotting
                    densityctrlfreq = 30,    # Frequency of density controls
                    maxiter = 10_000,        # When to terminate if no convergence
                    # Algorithmic parameters, defaults are normally good
                    Fscale = 1.2,            # Force function parameter (compression)
                    deltat = 0.2,            # Timestep
                    ttol = 0.1 * h0,         # Retriangulation tolerance
                    geps = 0.001 * h0,       # dfcn(p) tolerance for in/out points
                    dptol = geps,            # Termination tolerance
                    deps = sqrt(eps()) * h0  # Finite difference stepsize
                    )
    
    # Initializations
    pfix = unique(SF2.(pfix))
    nfix = length(pfix)
    pold = [SF2(Inf,Inf)]
    t, bars = SI3[], SI2[]
    converged = false
    
    # Initial nodes
    p = init_nodes(bbox, h0)
    p = nodes_rejection(p, pfix, dfcn, hfcn, geps)

    # Main loop
    for iter = 1:maxiter
        # If large relative movements, retriangulate and plot
        if maximum(norm.(p .- pold)) > ttol
            pold = copy(p)
            t, bars = retriangulate(p, dfcn, geps)
            plotting && display(simpplot(p, t))
        end

        # All bars and lengths
        barvec = barvectors(p, bars)
        L = norm.(barvec)
        L0 = desiredlengths(p, bars, L, hfcn, Fscale)
        
        # Density control - remove points that are too close to each other
        if iter % densityctrlfreq == 0
            p = density_control(p, L, L0, bars, nfix)
            pold = [SF2(Inf,Inf)]
            continue
        end

        # Compute forces and update nodes
        dp = deltat * total_node_forces(L, L0, barvec, bars, length(p), nfix)
        p .+= dp

        # Project outside points back to boundary
        d = project_nodes!(p, dfcn, deps)

        # Terminate if small (interior) node movements
        converged = maximum(norm.(dp[d.<-geps])) < dptol
        converged && break
    end
    
    converged || @warn "No convergence in maxiter=$maxiter iterations"
    plotting && display(simpplot(p, t))
    return DMesh(p, t)
end
