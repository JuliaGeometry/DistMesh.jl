"""
    distmesh
    3D Mesh Generator using Signed Distance Functions.
    Arguments:
        fdist:       Distance function
        fh:          Edge length function
        h:           Smallest edge length

    Returns:
        p:           Node positions
        t:           Triangle indices


    Example: Unit ball
        d(p) = sqrt(sum(p.^2))-1
        p,t = distmeshnd(d,huniform,0.2)
"""
function distmesh(fdist::Function,
                  fh::Union{Function,HUniform},
                  h::Number,
                  setup::AbstractDistMeshAlgorithm=DistMeshSetup();
                  origin=GeometryBasics.Point{3,Float64}(-1,-1,-1),
                  widths=GeometryBasics.Point{3,Float64}(2,2,2),
                  fix=nothing,
                  stats=false) where {VertType}
    # TODO: tetgen only handles Float64
    VT = GeometryBasics.Point{3,Float64}
    if isa(fix, Nothing)
        fp = nothing
    else
        fp = convert(Vector{VT}, fix)
    end
    o = VT(origin...)
    w = VT(widths...)
    distmesh(fdist, fh, h, setup, o, w, fp, Val(stats), VT)
end

"""
    DistMeshResult

A struct returned from the `distmesh` function that includes point, simplex,
and interation statistics.
"""
struct DistMeshResult{PT, TT, STATS}
    points::Vector{PT}
    tetrahedra::Vector{TT}
    stats::STATS
end

function distmesh(fdist::Function,
                  fh,
                  h::Number,
                  setup::DistMeshSetup,
                  origin,
                  widths,
                  fix,
                  ::Val{stats},
                  ::Type{VertType}) where {VertType, stats}

    geps=1e-1*h+setup.iso # parameter for filtering tets outside bounds and considering for max displacment of a node

    # static parameter info
    non_uniform = isa(fh, Function) # so we can elide fh calls
    have_fixed = !isa(fix, Nothing)

    #ptol=.001; ttol=.1; L0mult=1+.4/2^(dim-1); deltat=.2; geps=1e-1*h;

    # initialize Vertex Arrays
    # the first N points in the array correspond to'fix' points that do not move
    if have_fixed
        p = copy(fix)
    else
        p = VertType[]
    end

    pt_dists = eltype(VertType)[]

    # setup the initial point distribution specified in setup
    point_distribution!(fdist,p,pt_dists,h, setup, origin, widths, VertType)

    # Result struct for holding points, simplices, and iteration statistics
    result = DistMeshResult(p,
                            GeometryBasics.SimplexFace{4,Int32}[],
                            stats ? DistMeshStatistics() : nothing)

    # initialize arrays
    pair_set = Set{UInt64}()        # set used for ensure we have a unique set of edges
    pair = Tuple{Int32,Int32}[]                 # edge indices (Int32 since we use Tetgen)
    dp = zeros(VertType, length(p))             # displacement at each node
    bars = VertType[]                           # the vector of each edge
    L = eltype(VertType)[]                      # vector length of each edge
    L0 = non_uniform ? eltype(VertType)[] : nothing # desired edge length computed by dh (edge length function)
    maxmove = typemax(eltype(VertType))         # stores an iteration max movement for retriangulation

    # information on each iteration
    lcount = 0 # iteration counter
    triangulationcount = 0 # triangulation counter
    num_pairs = 0 # number of edge pairs

    @inbounds while true
        # if large move, retriangulation
        if maxmove>setup.ttol*h

            # compute a new delaunay triangulation
            retriangulate!(fdist, result, geps, setup, triangulationcount, pt_dists)

            num_pairs = tet_to_edges!(pair, pair_set, result.tetrahedra) # Describe each edge by a unique pair of nodes

            # resize arrays for new pair counts
            if triangulationcount == 0
                resize!(bars, length(result.tetrahedra)*6)
                resize!(L, length(result.tetrahedra)*6)
                non_uniform && resize!(L0, length(result.tetrahedra)*6)
            end

            triangulationcount += 1
            stats && push!(result.stats.retriangulations, lcount)
        end

        compute_displacements!(fh, dp, pair, num_pairs, L, L0, bars, result.points, setup, VertType)

        # Zero out forces on fix points
        if have_fixed
            for i in eachindex(fix)
                dp[i] = zero(VertType)
            end
        end

        # apply point forces and
        # bring outside points back to the boundary
        maxdp = typemin(eltype(VertType))
        maxmove = typemin(eltype(VertType))
        for i in eachindex(p)

            p0 = p[i] # store original point location
            p[i] = p[i].+setup.deltat.*dp[i] # apply displacements to points

            # Check if we are verifiably within the bounds and use this value
            # to avoid recomputing fdist. This increases performance greatly on
            # complex distance functions and large node counts
            move = sqrt(sum((p[i]-p0).^2))          # compute movement from the displacement
            d_est = pt_dists[i] + move              # apply the movement to our cache point
            d = d_est < -geps ? d_est : fdist(result.points[i]) # determine if we need correct or approximate distance
            pt_dists[i] = d                         # store distance

            if d < -geps
                maxdp = max(maxdp, setup.deltat*sqrt(sum(dp[i].^2)))
            end

            if d <= setup.iso
                maxmove = max(move,maxmove)
                continue
            end

            # bring points back to boundary if outside using central difference
            p[i] = p[i] .- centraldiff(fdist,p[i]).*(d+setup.iso)
            maxmove = max(sqrt(sum((p[i]-p0).^2)), maxmove)
            pt_dists[i] = setup.iso # ideally
        end

        # increment iteration counter
        lcount = lcount + 1

        # save iteration stats
        if stats
            push!(result.stats.maxmove,maxmove)
            push!(result.stats.maxdp,maxdp)
            min_v_edge, avg_v_edge, max_v_edge = volume_edge_stats(result.points,result.tetrahedra)
            push!(result.stats.min_volume_edge_ratio, min_v_edge)
            push!(result.stats.average_volume_edge_ratio, avg_v_edge)
            push!(result.stats.max_volume_edge_ratio, max_v_edge)
        end

        # Termination criterion
        if maxdp<setup.ptol*h
            return result
        end
    end
end

"""
    retriangulate!


Given a point set, generate a delaunay triangulation, and other requirements.
This includes:
    - Spatial sorting of points
    - Delaunay triangulation
    - Filtering of invalid tetrahedra outside the boundary
"""
function retriangulate!(fdist, result::DistMeshResult, geps, setup, triangulationcount, pt_dists)
    # use hilbert sort to improve cache locality of points
    if setup.sort && iszero(triangulationcount % setup.sort_interval)
        hilbertsort!(result.points, pt_dists)
    end

    t = result.tetrahedra
    p = result.points
    triangulation = delaunayn(p)
    t_d = triangulation.tetrahedra
    resize!(t, length(t_d))
    copyto!(t, t_d) # we need to copy since we have a shared reference with tetgen

    # average points to get mid point of each tetrahedra
    # if the mid point of the tetrahedra is outside of
    # the boundary we remove it.
    # TODO: this is an inlined filter call. Would be good to revert
    # TODO: can we use the point distance array to pass boundary points to
    #        tetgen so this call is no longer required?
    j = firstindex(t)
    @inbounds for ai in t
        t[j] = ai
        pm = (p[ai[1]].+p[ai[2]].+p[ai[3]].+p[ai[4]])./4
        j = ifelse(fdist(pm) <= -geps, nextind(t, j), j)
    end
    j <= lastindex(t) && resize!(t, j-1)
    nothing
end


function compute_displacements!(fh, dp, pair, num_pairs, L, L0, bars, p, setup,
    ::Type{VertType}) where {VertType}

    non_uniform = isa(typeof(L0), AbstractVector)

    # compute edge lengths (L) and adaptive edge lengths (L0)
    # Lp norm (p=3) is partially computed here
    Lsum = zero(eltype(L))
    L0sum = non_uniform ? zero(eltype(L0)) : length(pair)
    @inbounds for i in 1:num_pairs
        pb = pair[i]
        b1 = p[pb[1]]
        b2 = p[pb[2]]
        barvec = b1 - b2 # bar vector
        bars[i] = barvec
        L[i] = sqrt(sum(barvec.^2)) # length
        non_uniform && (L0[i] = fh((b1+b2)./2))
        Lsum = Lsum + L[i].^3
        non_uniform && (L0sum = L0sum + L0[i].^3)
    end

    # zero out force at each node
    @inbounds for i in eachindex(dp)
        dp[i] = zero(VertType)
    end

    # this is not hoisted correctly in the loop so we initialize here
    # finish computing the Lp norm (p=3)
    lscbrt = (1+(0.4/2^2))*cbrt(Lsum/L0sum)

    # Move mesh points based on edge lengths L and forces F
    @inbounds for i in 1:num_pairs
        if non_uniform && L[i] < L0[i]*lscbrt || L[i] < lscbrt
            L0_f = non_uniform ? L0[i].*lscbrt : lscbrt
            # compute force vectors
            F = setup.nonlinear ? (L[i]+L0_f)*(L0_f-L[i])/(2*L0_f) : L0_f-L[i]
            # edges are not allowed to pull, only repel
            FBar = bars[i].*F./L[i]
            # add the force vector to the node
            b1 = pair[i][1]
            b2 = pair[i][2]
            dp[b1] = dp[b1] .+ FBar
            dp[b2] = dp[b2] .- FBar
        end
    end
end
