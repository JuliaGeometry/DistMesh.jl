function distmeshnd(fdist,fh,h, ::Type{VertType}=GeometryBasics.Point{3,Float64}; origin=VertType(-1,-1,-1),
                                                                       widths=VertType(2,2,2), samples=_DEFAULT_SAMPLES) where {VertType}
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
    ptol=.001; ttol=.1; L0mult=1+.4/2^(dim-1); deltat=.1; geps=1e-1*h; deps=sqrt(eps())*h;

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

    nx, ny, nz = samples[1], samples[2], samples[3]


    # we subtract one from the length along each axis because
    # an NxNxN SDF has N-1 cells on each axis
    s = VertType(widths[1]/(nx-1), widths[2]/(ny-1), widths[3]/(nz-1))

    @inbounds for xi = 0:nx, yi = 0:ny, zi = 0:nz
        point = VertType(xi,yi,zi) .* s .+ origin
        fdist(point) <= geps && push!(p,point)
    end

    count=0;
    p0=fill(VertType(Inf),length(p))
    while true
        #% 3. Retriangulation by Delaunay
        # determine movements
        maxmove = -Inf
        tl = max(length(p),length(p0))
        for i in 1:tl
            maxmove = max(sqrt(sum((p[i]-p0[i]).^2)),maxmove)
        end
        if maxmove>ttol*h
            p0=copy(p)
            triangulation=delaunayn(p)
            t = copy(triangulation.tetrahedra)
            #pmid=zeros(size(t,1),dim);
            # average points to get mid point of each tetrahedra
            # TODO this is hardcoded for 3d
            #@show t[1:10], t[end-10:end]
            pmid = VertType[]
            for te in 1:length(t)
                #@show t[te]
                pm = sum(getindex(p,t[te]))/4
                push!(pmid,pm)
            end
            #pmid = [sum(getindex(p,te))/4 for te in eachindex(t)] # jl
            # for ii=1:dim+1
            #     pmid=pmid+p(t(:,ii),:)/(dim+1);
            # end
            #t=t(feval(fdist,pmid,varargin{:})<-geps,:);
            deletes = Int[]
            for i in eachindex(pmid)
                fdist(pmid[i]) < -geps && push!(deletes, i)
            end
            deleteat!(t, deletes)
            # % 4. Describe each edge by a unique pair of nodes
            pair=Vector{Tuple{Int,Int}}()
            # for ii=1:size(localpairs,1)
            #     pair=[pair;t(:,localpairs(ii,:))];
            # end
            for i in eachindex(t)
                for ep in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))
                    p1 = t[i][ep[1]]
                    p2 = t[i][ep[2]]
                    if p1 > p2
                        push!(pair, (p2,p1))
                    else
                        push!(pair, (p1,p2))
                    end
                end
            end
            unique!(pair)
            #pair=unique(sort(pair,2),'rows');
            #pair=munique(sort(pair,2)) #jl
            # % 5. Graphical output of the current mesh
            # if dim==2
            # trimesh(t,p(:,1),p(:,2),zeros(N,1))
            # view(2),axis equal,axis off,drawnow
            # elseif dim==3
            # if mod(count,5)==0
            #     simpplot(p,t,'p(:,2)>0');
            #     title(['Retriangulation #',int2str(count)])
            #     drawnow
            # end
            # else
            #     disp(sprintf('Retriangulation #%d',count))
            # end
            #count=count+1;
        end

        # 6. Move mesh points based on edge lengths L and forces F
        # bars=p(pair(:,1),:)-p(pair(:,2),:); # bar vector
        # L=sqrt(sum(bars.^2,2)); # length
        # L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
        # L0=L0*L0mult*(sum(L.^dim)/sum(L0.^dim))^(1/dim);
        # F=max(L0-L,0);
        # Fbar=[bars,-bars].*repmat(F./L,1,2*dim);
        # dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]),
        #                 ones(size(pair,1),1)*[1:dim,1:dim],
        #                 Fbar,N,dim));
        # dp(1:size(fix,1),:)=0;
        # p=p+deltat*dp;

        bars=[p[pb[1]]-p[pb[2]] for pb in pair] # bar vector
        L=[sqrt(sum(b.^2)) for b in bars] # length
        L0 = map(fh,[(p[pb[1]]+p[pb[2]])./2 for pb in pair])
        for i in 1:length(L0)
            L0[i] = L0mult*(sum(L[i].^dim)/sum(L0[i].^dim))^(1/dim)
        end
        F=[max(L0[i]-L[i],0) for i in eachindex(L0)]
        # TODO
        # Fbar=[bars,-bars].*repmat(F./L,1,2*dim)
        # dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]),
        #                 ones(size(pair,1),1)*[1:dim,1:dim],
        #                 Fbar,N,dim));
        # dp(1:size(fix,1),:)=0;
        # p=p+deltat*dp;

        #Fbar=[bars,-bars].*repmat(F./L,1,2*dim)
        FBar = F./L
        dp = fill(VertType(0), length(p))
        # sum up forces
        for i in eachindex(pair)
            b1 = pair[1]
            b2 = pair[2]
            p[b1] = p[b1] + FBar[i]
            p[b2] = p[b2] - FBar[i]
        end
        #dp(1:size(fix,1),:)=0;
        p=p.+deltat.*dp # apply displacements to points

        # 7. Bring outside points back to the boundary
        # d=feval(fdist,p,varargin{:}); ix=d>0;
        maxdp = -Inf
        for i in eachindex(p)
            d = fdist(p[i])
            d <= 0 && continue
            dx = (fdist(p[i].+VertType(deps,0,0)) + fdist(p[i].-VertType(deps,0,0)))/deps*2
            dy = (fdist(p[i].+VertType(0,deps,0)) + fdist(p[i].-VertType(0,deps,0)))/deps*2
            dz = (fdist(p[i].+VertType(0,0,deps)) + fdist(p[i].-VertType(0,0,deps)))/deps*2
            p[i] = p[i] - VertType(dz,dy,dz).*d
        end
        # gradd=zeros(sum(ix),dim);
        # for ii=1:dim
        #     a=zeros(1,dim);
        #     a(ii)=deps;
        #     d1x=feval(fdist,p(ix,:)+ones(sum(ix),1)*a,varargin{:});
        #     gradd(:,ii)=(d1x-d(ix))/deps;
        # end
        # p(ix,:)=p(ix,:)-d(ix)*ones(1,dim).*gradd;

        # 8. Termination criterion
        #maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));
        #if maxdp<ptol*h
        #     break
        #end
    end
end
