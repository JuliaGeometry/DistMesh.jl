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
function distmesh(fdist::Function,fh::Function,h::Number, setup::DistMeshSetup{T,ReTri}=DistMeshSetup(), ::Type{VertType}=GeometryBasics.Point{3,Float64}; origin=VertType(-1,-1,-1),
                                                                       widths=VertType(2,2,2),
                                                                       fix::Vector{VertType}=VertType[],
                                                                       stats=false,
                                                                       distribution=:regular) where {VertType, T, ReTri}

    ptol=setup.ptol
    L0mult=1+.4/2^2
    deltat=setup.deltat
    geps=1e-1*h+setup.iso
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

    if distribution == :regular
        simplecubic!(fdist, p, h, setup.iso, origin, widths, VertType)
    elseif distribution == :packed
        # face-centered cubic point distribution
        facecenteredcubic!(fdist, p, h, setup.iso, origin, widths, VertType)
    end
    dcount = 0
    lcount = 0

    # initialize arrays
    max_elt = typemax(eltype(VertType))
    pair = Tuple{Int32,Int32}[] # edge indices (Int32 since we use Tetgen)
    dp = fill(zero(VertType), length(p)) # force at each node
    bars = VertType[] # the vector of each edge
    L = eltype(VertType)[] # vector length of each edge
    L0 = eltype(VertType)[] # desired edge length computed by dh (edge length function)
    t = GeometryBasics.SimplexFace{4,Int32}[] # tetrahedra indices from delaunay triangulation
    maxmove = typemax(eltype(VertType)) # stores an iteration max movement for retriangulation

    # array for tracking quality metrics
    tris = NTuple{3,Int32}[] # array to store triangles used for quality checks
    qualities = eltype(VertType)[]
    #maxmoves = eltype(VertType)[]

    # information on each iteration
    statsdata = DistMeshStatistics()

    @inbounds while true
        # Retriangulation by Delaunay

        # if large move, retriangulation
        if ReTri <: RetriangulateMaxMove && maxmove>setup.retriangulation_criteria.ttol*h
            triangulation = delaunayn(p)
            t_d = triangulation.tetrahedra
            resize!(t, length(t_d))
            copyto!(t, t_d) # we need to copy since we have a shared reference with tetgen
            sort!(t) # sort tetrahedra so points are closer in mem

            # average points to get mid point of each tetrahedra
            # if the mid point of the tetrahedra is outside of
            # the boundary we remove it.
            if setup.droptets
                filter!(t) do i
                    pm = sum(getindex(p,i))/4
                    fdist(pm) <= -geps
                end
            end

            # 4. Describe each edge by a unique pair of nodes
            resize!(pair, length(t)*6)

            for i in eachindex(t)
                for ep in 1:6
                    p1 = t[i][tetpairs[ep][1]]
                    p2 = t[i][tetpairs[ep][2]]
                    if p1 > p2
                        pair[(i-1)*6+ep] = (p2,p1)
                    else
                        pair[(i-1)*6+ep] = (p1,p2)
                    end
                end
            end
            sort!(pair)
            unique!(pair)
            # resize arrays for new pair counts
            resize!(bars, length(pair))
            resize!(L, length(pair))
            resize!(L0, length(pair))

            stats && push!(statsdata.retriangulations, lcount)
        end

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
            deps = sqrt(eps(d))
            # use central difference
            dx = (fdist(p[i].+VertType(deps,0,0)) - fdist(p[i].-VertType(deps,0,0)))/(2deps)
            dy = (fdist(p[i].+VertType(0,deps,0)) - fdist(p[i].-VertType(0,deps,0)))/(2deps)
            dz = (fdist(p[i].+VertType(0,0,deps)) - fdist(p[i].-VertType(0,0,deps)))/(2deps)
            grad = VertType(dx,dy,dz) #normalize?
            # project back to boundary
            p[i] = p[i] - grad.*(d+setup.iso)
            maxmove = max(sqrt(sum((p[i]-p0).^2)),maxmove)
        end

        # increment iteration counter
        lcount = lcount + 1

        # save iteration stats
        if stats
            push!(statsdata.maxmove,maxmove)
            push!(statsdata.maxdp,maxdp)
            triangle_qualities!(tris,qualities,p,t)
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
