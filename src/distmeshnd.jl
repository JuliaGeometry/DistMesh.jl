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

function distmesh(fdist::Function,
                  fh,
                  h::Number,
                  setup::SetupTy,
                  origin,
                  widths,
                  fix,
                  ::Val{stats},
                  ::Type{VertType}) where {VertType, stats, SetupTy <: AbstractDistMeshAlgorithm}

    geps=1e-1*h+setup.iso # parameter for filtering tets outside bounds and considering for max displacment of a node

    # static parameter info
    non_uniform = isa(fh, Function) # so we can elide fh calls
    have_fixed = !isa(fix, Nothing)
    quality_mesh = isa(SetupTy, DistMeshQuality) # if we are using quality based retri and termination

    #ptol=.001; ttol=.1; L0mult=1+.4/2^(dim-1); deltat=.2; geps=1e-1*h;

    # initialize Vertex Arrays
    # the first N points in the array correspond to'fix' points that do not move
    if have_fixed
        p = copy(fix)
    else
        p = VertType[]
    end

    pt_dists = map(fdist, p) # cache to store point locations so we can minimize fdist calls

    # add points to p based on the initial distribution
    if setup.distribution === :regular
        simplecubic!(fdist, p, pt_dists, h, setup.iso, origin, widths, VertType)
    elseif setup.distribution === :packed
        # face-centered cubic point distribution
        facecenteredcubic!(fdist, p, pt_dists, h, setup.iso, origin, widths, VertType)
    end

    # initialize arrays
    pair_set = Set{Tuple{Int32,Int32}}()        # set used for ensure we have a unique set of edges
    pair = Tuple{Int32,Int32}[]                 # edge indices (Int32 since we use Tetgen)
    dp = zeros(VertType, length(p))             # displacement at each node
    bars = VertType[]                           # the vector of each edge
    L = eltype(VertType)[]                      # vector length of each edge
    L0 = eltype(VertType)[]                     # desired edge length computed by dh (edge length function)
    t = GeometryBasics.SimplexFace{4,Int32}[]   # tetrahedra indices from delaunay triangulation
    maxmove = typemax(eltype(VertType))         # stores an iteration max movement for retriangulation

    # arrays for tracking quality metrics
    if quality_mesh || stats
        tris = NTuple{3,Int32}[]        # triangles used for quality checks
        triset = Set{NTuple{3,Int32}}() # set for triangles to ensure uniqueness
        qualities = eltype(VertType)[]
        prev_mean_qual = typemax(eltype(VertType))
        new_mean_qual = typemin(eltype(VertType))
        min_qual = typemin(eltype(VertType))
    end

    if quality_mesh
        p_old = copy(p)
    end

    # information on each iteration
    statsdata = DistMeshStatistics()
    lcount = 0 # iteration counter
    triangulationcount = 0 # triangulation counter

    @inbounds while true
        # if large move, retriangulation
        # alternatively for quality metris we use the mean quality
        if (!quality_mesh && maxmove>setup.ttol*h) || quality_mesh && prev_mean_qual > new_mean_qual

            # if we are using quality meshing we revert the points
            # TODO: There might be places mean qual gets stuck
            if quality_mesh
                copyto!(p,p_old)
            end

            # use hilbert sort to improve cache locality of points
            # TODO: We can't sort if we have fixed
            if setup.sort && iszero(triangulationcount % setup.sort_interval)
                hilbertsort!(p)
            end

            # compute a new delaunay triangulation
            # we use the default random insertion method
            delaunayn!(fdist, p, t, geps, false)

            tet_to_edges!(pair, pair_set, t) # Describe each edge by a unique pair of nodes

            # resize arrays for new pair counts
            resize!(bars, length(pair))
            resize!(L, length(pair))
            non_uniform && resize!(L0, length(pair))

            # if the points were sorted we need to update the distance cache
            # TODO: We can't sort if we have fixed
            if setup.sort && iszero(triangulationcount % setup.sort_interval)
                for i in eachindex(p)
                    pt_dists[i] = fdist(p[i])
                end
            end

            # if we are using a quality based retriangulation, we need triangles
            if quality_mesh || stats
                # get triangles
                tets_to_tris!(tris, triset, t)
                resize!(qualities, length(tris))
            end

            stats && push!(statsdata.retriangulations, lcount)

            triangulationcount += 1
        end

        compute_displacements!(fh, dp, pair, L, L0, bars, p, setup, VertType)

        # Zero out forces on fix points
        if have_fixed
            for i in eachindex(fix)
                dp[i] = zero(VertType)
            end
        end

        # TODO: We don't really need displacmenets for quality mesh, but we can use it for stats
        maxdp = typemin(eltype(VertType))
        maxmove = typemin(eltype(VertType))
        # apply point forces and
        # bring outside points back to the boundary
        for i in eachindex(p)

            p0 = p[i] # store original point location
            p[i] = p[i].+setup.deltat.*dp[i] # apply displacements to points

            # Check if we are verifiably within the bounds and use this value
            # to avoid recomputing fdist. This increases performance greatly on
            # complex distance functions and large node counts
            move = sqrt(sum((p[i]-p0).^2))          # compute movement from the displacement
            d_est = pt_dists[i] + move              # apply the movement to our cache point
            d = d_est < -geps ? d_est : fdist(p[i]) # determine if we need correct or approximate distance
            pt_dists[i] = d                         # store distance

            # if we are not using quality improvements for retriangulation and termination,
            # track displacments and movements
            if d < -geps
                maxdp = max(maxdp, setup.deltat*sqrt(sum(dp[i].^2)))
            end

            if d <= setup.iso
                maxmove = max(move,maxmove)
                continue
            end

            # bring points back to boundary if outside using central difference
            p[i] = p[i] .- centraldiff(fdist,p[i]).*(d+setup.iso)
            !quality_mesh && (maxmove = max(sqrt(sum((p[i]-p0).^2)), maxmove))
            pt_dists[i] = setup.iso # ideally
        end

        # increment iteration counter
        lcount = lcount + 1

        # some setup for quality mesh or stats collection
        if quality_mesh || stats
            triangle_qualities!(tris, qualities, p, t)
            min_qual = minimum(qualities)
            mean_qual = sum(qualities)/length(qualities)
            prev_mean_qual = new_mean_qual
            new_mean_qual = mean_qual
        end
        # save iteration stats
        if stats
            push!(statsdata.maxmove,maxmove)
            push!(statsdata.maxdp,maxdp)
            sort!(qualities) # sort for median calc and robust summation
            max_qual = maximum(qualities)
            push!(statsdata.average_qual, mean_qual)
            # TODO: We can use Statistic Median! here without doing a full sort
            push!(statsdata.median_qual, qualities[round(Int,length(qualities)/2)])
            push!(statsdata.minimum_qual, min_qual)
            push!(statsdata.maximum_qual, max_qual)
        end

        # Termination criterion
        if (!quality_mesh && maxdp<setup.ptol*h) || (quality_mesh && min_qual < setup.min_quality)
            return p, t, statsdata
        end
    end
end
