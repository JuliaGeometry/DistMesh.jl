function mkt2t(t)
    #MKT2T  Compute element connectivities from element indices.
    #   [T2T,T2N]=MKT2T(T)

    #   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

    nt=size(t,1)
    dim=size(t,2)-1

    # switch dim
    #  case 1
    #   edges=[t(:,2)
    #          t(:,1)];
    #  case 2
    #   edges=[t(:,[2,3])
    #          t(:,[3,1])
    #          t(:,[1,2])];
    #  case 3
    #   edges=[t(:,[2,3,4])
    #          t(:,[3,4,1])
    #          t(:,[4,1,2])
    #          t(:,[1,2,3])];
    # end

    if dim == 1
        edges=[t[:,2]
                t[:,1]];
    elseif dim == 2
        edges=[t[:,[2,3]]
                t[:,[3,1]]
                t[:,[1,2]]];
    elseif dim == 3
        edges=[t[:,[2,3,4]]
                t[:,[3,4,1]]
                t[:,[4,1,2]]
                t[:,[1,2,3]]];
    end #ok

    # ts=[repmat(int32(1:nt),1,dim+1); kron(int32(1:(dim+1)),ones(1,nt,'int32'))]';
    # this might be done with just "dim" number of generators or something.
    # we want 1:nt (cross) 1:dim+1
    ts=hcat(repeat(1:nt,dim+1), kron(1:(dim+1),ones(Int, 1,nt)')) # ok

    # edges=sort(edges,2);
    edges=sort(edges,dims=2) # ok
    jx=munique(edges)
    @show length(jx), jx[1:15]

    # [jx,ix]=sort(jx);
    ix = sortperm(jx)
    jx = sort(jx)
    @show length(jx), jx[1:15,:] # okay
    @show length(ix), ix[1:15,:] # okay

    # ts=ts(ix,:);
    ts=ts[ix,:]

    # ix=find(diff(jx)==0);
    djx = diff(jx)
    @show typeof(jx), typeof(djx)
    ix=findall(iszero, djx) #okay
    # ts1=ts(ix,:);
    @show length(ix), ix[1:15,:] # okay
    ts1=ts[ix,:]
    @show length(ts1), ts1[1:15,:]
    # ts2=ts(ix+1,:);
    ts2=ts[ix.+1,:]
    @show length(ts2), ts2[1:15,:]

    # t2t=zeros(nt,dim+1,'int32');
    # t2t(ts1(:,1)+nt*(ts1(:,2)-1))=ts2(:,1);
    # t2t(ts2(:,1)+nt*(ts2(:,2)-1))=ts1(:,1);
    t2t=zeros(Int, nt,dim+1)
    @show typeof(ts1), typeof(ts2)
    t2t[ts1[:,1]+nt*(ts1[:,2].-1)]=ts2[:,1]
    t2t[ts2[:,1]+nt*(ts2[:,2].-1)]=ts1[:,1]
    @show length(t2t), t2t[1:15,:]

    # if nargout>=2
    #   t2n=zeros(nt,dim+1,'int32');
    #   t2n(ts1(:,1)+nt*(ts1(:,2)-1))=ts2(:,2);
    #   t2n(ts2(:,1)+nt*(ts2(:,2)-1))=ts1(:,2);
    # end

    t2n=zeros(Int, nt,dim+1);
    t2n[ts1[:,1]+nt*(ts1[:,2].-1)]=ts2[:,2]
    t2n[ts2[:,1]+nt*(ts2[:,2].-1)]=ts1[:,2]
    t2t, t2n
end