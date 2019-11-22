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
                  setup::DistMeshSetup=DistMeshSetup();
                  origin=GeometryBasics.Point{3,Float64}(-1,-1,-1),
                  widths=GeometryBasics.Point{3,Float64}(2,2,2),
                  fix=nothing,
                  stats=false) where {VertType}
    if isa(fix, Nothing)
        VT = promote_type(typeof(origin), typeof(widths))
        fp = nothing
    else
        VT = promote_type(typeof(origin), typeof(widths), eltype(fix))
        fp = convert(Vector{VT}, fix)
    end
    o = VT(origin...)
    w = VT(widths...)
    distmesh(fdist, fh, h, setup, o, w, fp, Val(stats), VT)
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

    geps=1e-1*h+setup.iso

    # static parameter info
    non_uniform = isa(fh, Function) # so we can elide fh calls
    have_fixed = !isa(fix, Nothing)

    deps = sqrt(eps(eltype(VertType))) # epsilon for computing central difference
    #ptol=.001; ttol=.1; L0mult=1+.4/2^(dim-1); deltat=.2; geps=1e-1*h;

    # # % 2. Remove points outside the region, apply the rejection method
    # p=p(feval(fdist,p,varargin{:})<geps,:);
    # r0=feval(fh,p);
    # p=[fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];
    # N=size(p,1);

    # initialize Vertex Arrays
    # the first N points in the array correspond to
    # 'fix' points that do not move
    if have_fixed
        p = copy(fix)
    else
        p = VertType[]
    end

    # add points to p based on the initial distribution
    if setup.distribution === :regular
        simplecubic!(fdist, p, h, setup.iso, origin, widths, VertType)
    elseif setup.distribution === :packed
        # face-centered cubic point distribution
        facecenteredcubic!(fdist, p, h, setup.iso, origin, widths, VertType)
    end

    pt_dists = map(fdist, p) # cache to store point locations so we can minimize fdist calls

    # initialize arrays
    pair_set = Set{Tuple{Int32,Int32}}() # set used for ensure we have a unique set of edges
    pair = Tuple{Int32,Int32}[] # edge indices (Int32 since we use Tetgen)
    dp = fill(zero(VertType), length(p)) # force at each node
    bars = VertType[] # the vector of each edge
    L = eltype(VertType)[] # vector length of each edge
    non_uniform && (L0 = eltype(VertType)[]) # desired edge length computed by dh (edge length function)
    t = GeometryBasics.SimplexFace{4,Int32}[] # tetrahedra indices from delaunay triangulation
    maxmove = typemax(eltype(VertType)) # stores an iteration max movement for retriangulation

    # array for tracking quality metrics
    tris = NTuple{3,Int32}[] # array to store triangles used for quality checks
    triset = Set{NTuple{3,Int32}}() # set for triangles to ensure uniqueness
    qualities = eltype(VertType)[]
    #maxmoves = eltype(VertType)[]

    # information on each iteration
    statsdata = DistMeshStatistics()
    lcount = 0 # iteration counter

    @inbounds while true
        # Retriangulation by Delaunay

        # if large move, retriangulation
        if maxmove>setup.ttol*h
            # compute a new delaunay triangulation
            delaunayn!(fdist, p, t, geps)

            # 4. Describe each edge by a unique pair of nodes
            tet_to_edges!(pair, pair_set, t)

            # resize arrays for new pair counts
            resize!(bars, length(pair))
            resize!(L, length(pair))
            non_uniform && resize!(L0, length(pair))

            stats && push!(statsdata.retriangulations, lcount)
        end

        # 6. Move mesh points based on edge lengths L and forces F
        Lsum = zero(eltype(L))
        L0sum = non_uniform ? zero(eltype(L0)) : length(pair)
        for i in eachindex(pair)
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
        for i in eachindex(dp)
            dp[i] = zero(VertType)
        end

        # this is not hoisted correctly in the loop
        # so we initialize here
        # this is the Lp norm (p=3)
        lscbrt = (1+(0.4/2^2))*cbrt(Lsum/L0sum)

        for i in eachindex(pair)
            L0_f = non_uniform ? L0[i].*lscbrt : lscbrt
            # compute force vectors
            F = L0_f-L[i]
            # edges are not allowed to pull, only repel
            if F > zero(eltype(L))
                FBar = bars[i].*F./L[i]
                # add the force vector to the node
                b1 = pair[i][1]
                b2 = pair[i][2]
                dp[b1] = dp[b1] .+ FBar
                dp[b2] = dp[b2] .- FBar
            end
        end

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

            # check if we are verifiably within the bounds and use this value
            # to avoid recomputing fdist
            move = sqrt(sum((p[i]-p0).^2))
            d_est = pt_dists[i] + move # apply the movement to our cache point
            d = d_est < -geps ? d_est : fdist(p[i])
            pt_dists[i] = d # store the correct or approximate distance from the function

            if d < -geps
                maxdp = max(maxdp, setup.deltat*sqrt(sum(dp[i].^2)))
            end

            if d <= setup.iso
                maxmove = max(move,maxmove) # determine movements
                continue
            end

            # bring points back to boundary if outside
            #deps = sqrt(eps(d)) # this is quite an expensive call, we use a constant initialized in the beginning
            # use central difference and project back to boundary
            p[i] = p[i] .- centraldiff(fdist,p[i]).*(d+setup.iso)
            maxmove = max(sqrt(sum((p[i]-p0).^2)), maxmove)
            pt_dists[i] = setup.iso # ideally
        end

        # increment iteration counter
        lcount = lcount + 1

        # save iteration stats
        if stats
            push!(statsdata.maxmove,maxmove)
            push!(statsdata.maxdp,maxdp)
            triangle_qualities!(tris,triset,qualities,p,t)
            sort!(qualities) # sort for median calc and robust summation
            push!(statsdata.average_qual, sum(qualities)/length(qualities))
            push!(statsdata.median_qual, qualities[round(Int,length(qualities)/2)])
        end

        # 8. Termination criterion
        if maxdp<setup.ptol*h
            return p, t, statsdata
        end
    end
end
