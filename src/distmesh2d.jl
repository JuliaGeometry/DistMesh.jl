########################################################################
# Internal Type Aliases

const Point2d = SVector{2, Float64}
const Index2 = SVector{2, Int32}   # For Edges
const Index3 = SVector{3, Int32}   # For Triangles

########################################################################
# Utility functions

barvectors(p, bars) = [ p[bar[1]] - p[bar[2]] for bar in bars ]

function init_nodes(bbox, h0)
    # Create initial distribution in bounding box (equilateral triangles)
    xx = bbox[1][1]:h0:bbox[2][1]
    yy = bbox[1][2]:h0*√3/2:bbox[2][2]
    p = [ Point2d(x + iseven(iy) * h0/2, y) for x in xx for (iy,y) in enumerate(yy) ]
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
    bars = [ tt[ix] for tt in t for ix in Index2[(1,2),(2,3),(3,1)] ]
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
    Ftot = [ Point2d(0.0,0.0) for ip = 1:np ]
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
    numgrad(f,x,fx) = (Point2d(f(x + Point2d(deps,0)), f(x + Point2d(0,deps))) .- fx) / deps
    dgrad = [ numgrad(dfcn, p[i], d[i]) for i in ix ]
    @. p[ix] -= d[ix] * dgrad / norm(dgrad)^2
    d
end

########################################################################
# Main function

"""
    distmesh2d(dfcn, hfcn, h0, bbox, pfix=[]; kwargs...) -> (p, t)

Generate a 2D unstructured triangular mesh using a signed distance function.

This function implements the DistMesh algorithm (Persson/Strang), which treats the mesh generation 
as a physical equilibrium problem. A system of truss bars (edges) is relaxed until the force 
equilibrium is reached, constrained by the signed distance function `dfcn` to stay within the domain.

# Arguments
- `dfcn::Function`: Signed distance function `d(p)`. Returns negative values inside the region, 
  positive outside. Input `p` is an `(N,2)` array.
- `hfcn::Function`: Element size function `h(p)`. Returns target edge length at point `p`.
- `h0::Real`: Initial nominal edge length (scaling factor for `hfcn`).
- `bbox`: Bounding box tuple `((xmin, ymin), (xmax, ymax))` defining the initial grid generation area.
- `pfix::Vector`: (Optional) List of fixed node positions that must be part of the mesh (e.g., corners).

# Keywords
- `plotting::Bool = false`: Enable live visualization of the relaxation process (requires GLMakie).
- `maxiter::Int = 1000`: Maximum number of relaxation iterations.
- `ttol::Real`: Tolerance for retriangulation (defaults to `0.1 * h0`).
- `dptol::Real`: Termination tolerance for node movement (defaults to `0.001 * h0`).

# Returns
- `p::Vector{Point2d}`: The node positions.
- `t::Vector{Index3}`: The triangle connectivity indices.

# Examples

**Uniform Mesh on a Unit Circle**
```julia
using DistMesh
# Distance function for a circle of radius 1
fd = p -> sqrt.(sum(p.^2, dims=2)) .- 1.0
# Uniform mesh size
fh = p -> ones(size(p,1))

p, t = distmesh2d(fd, fh, 0.2, ((-1,-1), (1,1)), plotting=false)
```

**Rectangle with Circular Hole (Refined at Boundary)**
```julia
using DistMesh

# Define distance function: Rectangle minus Circle
function fd_hole(p)
    d_rect = drectangle(p, -1, 1, -1, 1)
    d_circ = dcircle(p, 0, 0, 0.5)
    return ddiff(d_rect, d_circ)
end

# Size function: Refine near the hole
fh_hole(p) = 0.05 .+ 0.3 * dcircle(p, 0, 0, 0.5)

bbox = ((-1,-1), (1,1))
pfix = [(-1,-1), (-1,1), (1,-1), (1,1)] # Fix corners

p, t = distmesh2d(fd_hole, fh_hole, 0.05, bbox, pfix)
```
"""
function distmesh2d(dfcn, hfcn, h0, bbox, pfix=Point2d[];
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
    pfix = unique(Point2d.(pfix))
    nfix = length(pfix)
    pold = [Point2d(Inf,Inf)]
    t, bars = Index3[], Index2[]
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
            plotting && live_plot(DMesh(p,t))
        end

        # All bars and lengths
        barvec = barvectors(p, bars)
        L = norm.(barvec)
        L0 = desiredlengths(p, bars, L, hfcn, Fscale)
        
        # Density control - remove points that are too close to each other
        if iter % densityctrlfreq == 0
            p = density_control(p, L, L0, bars, nfix)
            pold = [Point2d(Inf,Inf)]
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

    msh = DMesh(p, t)
    plotting && live_plot(msh)
    return msh
end
