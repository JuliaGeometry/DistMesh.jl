function distmeshnd(fdist,fh,h, ::Type{VertType}=GeometryBasics.Point{3,Float64}; origin=VertType(-1,-1,-1),
                                                                       widths=VertType(2,2,2)) where {VertType}
    # %DISTMESHND N-D Mesh Generator using Distance Functions.
    # %   [P,T]=DISTMESHND(FDIST,FH,H,BOX,FIX,FDISTPARAMS)
    # %
    # %      P:           Node positions (NxNDIM)
    # %      T:           Triangle indices (NTx(NDIM+1))
    # %      FDIST:       Distance function
    # %      FH:          Edge length function
    # %      H:           Smallest edge length
    # %      BOX:         Bounding box [xmin,xmax;ymin,ymax; ...] (NDIMx2)
    # %      FIX:         Fixed node positions (NFIXxNDIM)
    # %      FDISTPARAMS: Additional parameters passed to FDIST
    # %
    # %   Example: Unit ball
    # %      dim=3;
    # %      d=inline('sqrt(sum(p.^2,2))-1','p');
    # %      [p,t]=distmeshnd(d,@huniform,0.2,[-ones(1,dim);ones(1,dim)],[]);
    # %
    # %   See also: DISTMESH2D, DELAUNAYN, TRIMESH, MESHDEMOND.

    # %   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

    dim=length(VertType)
    @show dim
    ptol=.001; ttol=.1; L0mult=1+.4/2^(dim-1); deltat=.2; geps=1e-1*h;

    # # 1. Create initial distribution in bounding box
    # if dim==1
    #     p=(box(1):h:box(2))';
    # else
    #     cbox=cell(1,dim);
    #     for ii=1:dim
    #         cbox{ii}=box(1,ii):h:box(2,ii);
    #     end
    #     pp=cell(1,dim);
    #     [pp{:}]=ndgrid(cbox{:});
    #     p=zeros(prod(size(pp{1})),dim);
    #     for ii=1:dim
    #         p(:,ii)=pp{ii}(:);
    #     end
    # end

    # # % 2. Remove points outside the region, apply the rejection method
    # p=p(feval(fdist,p,varargin{:})<geps,:);
    # r0=feval(fh,p);
    # p=[fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];
    # N=size(p,1);

    p = VertType[]

    # TODO try to bo back to uniform distribution here
    samples = round(reduce(*,widths./(h*1.5)))
    # for _ in 1:samples
    #     point = rand(VertType).*widths + origin
    #     @show point
    #     fdist(point) < -h && push!(p,point)
    # end
    # we subtract one from the length along each axis because
    # an NxNxN SDF has N-1 cells on each axis

    @inbounds for xi = origin[1]:h:(origin[1]+widths[1]), yi = origin[2]:h:(origin[2]+widths[2]), zi = origin[3]:h:(origin[3]+widths[3])
        point = VertType(xi,yi,zi)
        fdist(point) < 0 && push!(p,point)
    end
    dcount = 0
    lcount = 0

    # initialize arrays
    p0=fill(VertType(Inf),length(p))
    pair = Tuple{Int,Int}[]
    dp = fill(VertType(0), length(p))
    bars = VertType[]
    L = eltype(VertType)[]
    L0 = eltype(VertType)[]
    F = eltype(VertType)[]
    FBar = VertType[]
    t = GeometryBasics.SimplexFace{4,Int32}[]

    # makie viz
    ls = Pair{VertType,VertType}[]

    while true
        #@show dcount, lcount
        #% 3. Retriangulation by Delaunay
        # determine movements
        maxmove = -Inf
        tl = max(length(p),length(p0))
        for i in 1:tl
            maxmove = max(sqrt(sum((p[i]-p0[i]).^2)),maxmove)
        end
        #@show maxmove, ttol*h
        if maxmove>ttol*h
            triangulation=delaunayn(p)
            t = copy(triangulation.tetrahedra)
            # average points to get mid point of each tetrahedra
            # if the mid point of the tetrahedra is outside of
            # the boundary we remove it.
            # TODO this is hardcoded for 3d
            deletes = Int[]
            for i in eachindex(t)
                pm = sum(getindex(p,t[i]))/4
                fdist(pm) > -geps && push!(deletes, i)
            end
            deleteat!(t, deletes)
            # 4. Describe each edge by a unique pair of nodes
            pair=Vector{Tuple{Int,Int}}(undef, length(t)*6)

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
            unique!(pair)
            # resize arrays for new pair counts
            resize!(bars, length(pair))
            resize!(L, length(pair))
            resize!(L0, length(pair))
            resize!(F, length(pair))
            resize!(FBar, length(pair))

            # makie vis
            if dcount%5 == 0
                resize!(ls, length(pair))
                for i = 1:length(pair)
                    ls[i] = p[pair[i][1]] => p[pair[i][2]]
                end
                scene = Makie.linesegments(ls)
                display(scene)
                sleep(0.01)
            end
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
        L0 = L0.*L0mult.*(Lsum/L0sum)^(1/dim)
        for i in eachindex(L0)
            F[i] = max(L0[i]-L[i],zero(eltype(L0)))
        end

        # compute bar forces
        for i in eachindex(L0)
            FBar[i] = bars[i].*F[i]./L[i]
        end
        # zero out force at each node
        for i in eachindex(dp)
            dp[i] = zero(VertType)
        end
        # sum up forces
        for i in eachindex(pair)
            b1 = pair[i][1]
            b2 = pair[i][2]
            dp[b1] = dp[b1] .+ FBar[i]
            dp[b2] = dp[b2] .- FBar[i]
        end
        p0=copy(p)
        # TODO apply fixed points
        #dp(1:size(fix,1),:)=0;
        p=p.+deltat.*dp # apply displacements to points
        # 7. Bring outside points back to the boundary
        maxdp = -Inf
        for i in eachindex(p)
            d = fdist(p[i])
            if d < -geps
                maxdp = max(maxdp, deltat*sqrt(sum(dp[i].^2)))
            end
            d <= 0 && continue
            deps = sqrt(eps(d))
            # use central difference
            dx = (fdist(p[i].+VertType(deps,0,0)) - fdist(p[i].-VertType(deps,0,0)))/(2deps)
            dy = (fdist(p[i].+VertType(0,deps,0)) - fdist(p[i].-VertType(0,deps,0)))/(2deps)
            dz = (fdist(p[i].+VertType(0,0,deps)) - fdist(p[i].-VertType(0,0,deps)))/(2deps)
            grad = VertType(dx,dy,dz) #normalize?
            # project back to boundary
            p[i] = p[i] - grad.*d
        end
        lcount = lcount + 1
        # 8. Termination criterion
        #@show maxdp, ptol*h
        if maxdp<ptol*h
            @show sum(L)/length(L)
            return p, t
        end
    end
end
