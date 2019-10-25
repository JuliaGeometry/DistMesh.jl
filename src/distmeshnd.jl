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
function distmesh(fdist::Function,fh::Function,h::Number, ::Type{VertType}=GeometryBasics.Point{3,Float64}; origin=VertType(-1,-1,-1),
                                                                       widths=VertType(2,2,2),
                                                                       fix::Vector{VertType}=VertType[],
                                                                       vis=true,
                                                                       distribution=:regular) where {VertType}

    dim=length(VertType)
    ptol=.001; ttol=0.02; L0mult=1+.4/2^(dim-1); deltat=0.05; geps=1e-1*h;
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
        @inbounds for xi = origin[1]:h:(origin[1]+widths[1]), yi = origin[2]:h:(origin[2]+widths[2]), zi = origin[3]:h:(origin[3]+widths[3])
            point = VertType(xi,yi,zi)
            fdist(point) < 0 && push!(p,point)
        end
    elseif distribution == :packed
        # face-centered cubic point distribution
        r = h/2
        counts = round.(widths./h).+2
        @inbounds for xi = -1:Int(counts[1]), yi = -1:Int(counts[2]), zi = -1:Int(counts[3])
            point = VertType(2xi+((yi+zi)%2), sqrt(3)*(yi+(zi%2)/3),2*sqrt(6)*zi/3).*r + origin
            fdist(point) < 0 && push!(p,point)
        end
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

    # makie viz
    #ls = Pair{VertType,VertType}[]

    @inbounds while true
        # Retriangulation by Delaunay

        # if large move, retriangulation
        if maxmove>ttol*h
            triangulation=delaunayn(p)
            t_d = triangulation.tetrahedra
            resize!(t, length(t_d))
            for i in eachindex(t_d)
                t[i] = t_d[i]
            end
            # average points to get mid point of each tetrahedra
            # if the mid point of the tetrahedra is outside of
            # the boundary we remove it.
            # TODO this is hardcoded for 3d
            filter!(t) do i
                pm = sum(getindex(p,i))/4
                fdist(pm) <= -geps
            end

            # 4. Describe each edge by a unique pair of nodes
            pair=resize!(pair, length(t)*6)

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

            # makie vis
            # if vis
            #     if dcount%5 == 0
            #         resize!(ls, length(pair))
            #         for i = 1:length(pair)
            #             ls[i] = p[pair[i][1]] => p[pair[i][2]]
            #         end
            #         scene = Makie.linesegments(ls)
            #         display(scene)
            #         sleep(0.01)
            #     end
            # end
            dcount=dcount+1
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
            Lsum = Lsum + L[i].^dim
            L0sum = L0sum + L0[i].^dim
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

            if d <= 0
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
            p[i] = p[i] - grad.*d
            maxmove = max(sqrt(sum((p[i]-p0).^2)),maxmove)
        end
        lcount = lcount + 1
        # 8. Termination criterion
        #@show maxdp, ptol*h
        if maxdp<ptol*h
            @show dcount, lcount
            @show sum(L)/length(L)
            return p, t
        end
    end
end
