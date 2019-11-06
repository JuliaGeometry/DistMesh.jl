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
function distmesh(fdist::Function,fh::Function,h::Number, setup::DistMeshSetup{T}=DistMeshSetup(), ::Type{VertType}=GeometryBasics.Point{3,Float64}; origin=VertType(-1,-1,-1),
                                                                       widths=VertType(2,2,2),
                                                                       fix::Vector{VertType}=VertType[],
                                                                       stats=false) where {VertType, T}

    ptol=setup.ptol
    L0mult=1+.4/2^2
    deltat=setup.deltat
    geps=1e-1*h+setup.iso

    deps = sqrt(eps(T)) # epsilon for computing central difference
    #ptol=.001; ttol=.1; L0mult=1+.4/2^(dim-1); deltat=.2; geps=1e-1*h;

    # # % 2. Remove points outside the region, apply the rejection method
    # p=p(feval(fdist,p,varargin{:})<geps,:);
    # r0=feval(fh,p);
    # p=[fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];
    # N=size(p,1);

    # initialize Vertex Arrays
    # the first N points in the array correspond to
    # 'fix' points that do not move
    p = copy(fix)

    if setup.initial_points === :regular
        simplecubic!(fdist, p, h, setup.iso, origin, widths, VertType)
    elseif setup.intial_points === :packed
        # face-centered cubic point distribution
        facecenteredcubic!(fdist, p, h, setup.iso, origin, widths, VertType)
    end
    lcount = 0

    # initialize arrays
    pair_set = Set{Tuple{Int32,Int32}}() # set used for ensure we have a unique set of edges
    pair = Tuple{Int32,Int32}[] # edge indices (Int32 since we use Tetgen)
    dp = fill(zero(VertType), length(p)) # force at each node
    bars = VertType[] # the vector of each edge
    L = eltype(VertType)[] # vector length of each edge
    L0 = eltype(VertType)[] # desired edge length computed by dh (edge length function)
    t = GeometryBasics.SimplexFace{4,Int32}[] # tetrahedra indices from delaunay triangulation
    maxmove = typemax(eltype(VertType)) # stores an iteration max movement for retriangulation

    # array for tracking quality metrics
    tris = NTuple{3,Int32}[] # array to store triangles used for quality checks
    triset = Set{NTuple{3,Int32}}() # set for triangles to ensure uniqueness
    qualities = eltype(VertType)[]
    maxmoves = eltype(VertType)[]

    # information on each iteration
    statsdata = DistMeshStatistics()
    last_retri = 0 # iterations since last retriangulation

    @inbounds while true
        # Retriangulation by Delaunay
        if length(maxmoves) == setup.maxmove_delta
            popfirst!(maxmoves)
        end
        lcount > 0 && push!(maxmoves, maxmove) # just dont include the very first iteration
        # if large move, retriangulation
        if lcount < 7 && maxmove > setup.ttol*h || last_retri > setup.maxmove_delta_delay && maxmove > sum(maxmoves)/length(maxmoves)
            triangulation = delaunayn(p)
            t_d = triangulation.tetrahedra
            resize!(t, length(t_d))
            copyto!(t, t_d) # we need to copy since we have a shared reference with tetgen

            # average points to get mid point of each tetrahedra
            # if the mid point of the tetrahedra is outside of
            # the boundary we remove it.
            if setup.droptets
                j = firstindex(t)
                for ai in t
                    t[j] = ai
                    pm = (p[ai[1]].+p[ai[2]].+p[ai[3]].+p[ai[4]])./4
                    j = ifelse(fdist(pm) <= -geps, nextind(t, j), j)
                end
                j <= lastindex(t) && resize!(t, j-1)
            end

            # 4. Describe each edge by a unique pair of nodes
            empty!(pair_set)
            for i in eachindex(t)
                for ep in 1:6
                    p1 = t[i][tetpairs[ep][1]]
                    p2 = t[i][tetpairs[ep][2]]
                    push!(pair_set, p1 > p2 ? (p2,p1) : (p1,p2))
                end
            end
            resize!(pair, length(pair_set))
            # copy pair set to array since sets are not sortable
            i = 1
            for elt in pair_set
                pair[i] = elt
                i = i + 1
            end

            # sort the edge pairs for better point lookup
            sort!(pair)
            # resize arrays for new pair counts
            resize!(bars, length(pair))
            resize!(L, length(pair))
            resize!(L0, length(pair))

            empty!(maxmoves)
            last_retri = 0
            stats && push!(statsdata.retriangulations, lcount)
        end

        last_retri = last_retri + 1

        # 6. Move mesh points based on edge lengths L and forces F
        Lsum = zero(eltype(L))
        L0sum = zero(eltype(L0))
        for i in eachindex(pair)
            pb = pair[i]
            b1 = p[pb[1]]
            b2 = p[pb[2]]
            barvec = b1 - b2 # bar vector
            bars[i] = barvec
            L[i] = sqrt(sum(barvec.^2)) # length
            L0[i] = fh((b1+b2)./2)
            Lsum = Lsum + L[i].^3
            L0sum = L0sum + L0[i].^3
        end

        # zero out force at each node
        for i in eachindex(dp)
            dp[i] = zero(VertType)
        end

        for i in eachindex(pair)
            L0_f = L0[i].*L0mult.*cbrt(Lsum/L0sum)
            # compute force vectors
            F = max(L0_f-L[i],zero(eltype(L0)))
            FBar = bars[i].*F./L[i]
            # add the force vector to the node
            b1 = pair[i][1]
            b2 = pair[i][2]
            dp[b1] = dp[b1] .+ FBar
            dp[b2] = dp[b2] .- FBar
        end

        # Zero out forces on fix points
        for i in eachindex(fix)
            dp[i] = zero(VertType)
        end

        # apply point forces and
        # bring outside points back to the boundary
        maxdp = typemin(eltype(VertType))
        maxmove = typemin(eltype(VertType))
        for i in eachindex(p)

            p0 = p[i] # store original point location
            p[i] = p[i].+deltat.*dp[i] # apply displacements to points

            d = fdist(p[i])

            if d < -geps
                maxdp = max(maxdp, deltat*sqrt(sum(dp[i].^2)))
            end

            if d <= setup.iso
                maxmove = max(sqrt(sum((p[i]-p0).^2)),maxmove) # determine movements
                continue
            end

            # bring points back to boundary if outside
            #deps = sqrt(eps(d)) # this is quite an expensive call, we use a constant initialized in the beginning
            # use central difference
            dx = (fdist(p[i].+VertType(deps,0,0)) - fdist(p[i].-VertType(deps,0,0)))/(2deps)
            dy = (fdist(p[i].+VertType(0,deps,0)) - fdist(p[i].-VertType(0,deps,0)))/(2deps)
            dz = (fdist(p[i].+VertType(0,0,deps)) - fdist(p[i].-VertType(0,0,deps)))/(2deps)
            grad = VertType(dx,dy,dz) #normalize?
            # project back to boundary
            p[i] = p[i] .- grad.*(d+setup.iso)
            maxmove = max(sqrt(sum((p[i]-p0).^2)),maxmove)
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
        if maxdp<ptol*h
            return p, t, statsdata
        end
    end
end
