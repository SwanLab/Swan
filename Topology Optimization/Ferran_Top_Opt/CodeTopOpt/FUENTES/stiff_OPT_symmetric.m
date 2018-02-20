<<<<<<< HEAD
function [StifMat,post] = stiff_OPT_symmetric(dim,element,problembsc,coordinatesn,coordinatesa,gamma)
nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

% lnods = zeros(nnode,nelem);
idx   = zeros(nnode*nunkn,nelem);
ptype = problembsc.problemtype;
ftype = problembsc.phisical_type;
StifMat = sparse(nndof,nndof);
% StifMat_gauss = repmat({sparse(nndof,nndof)},1,element.ngaus);

lnods = element.conectivities';
% for i=1:nnode
%     lnods(i,:)= element.conectivities(:,i);
% end

etype = element.type;
for inode=1:nnode
    for idime=1:nunkn
        idx(nunkn*inode-nunkn+idime,:) = nunkn.*lnods(inode,:)-nunkn+idime;
    end
end
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
post.chi = zeros(ngaus,nelem);

% tic
for igaus=1:ngaus
     [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [Bmat] = cal_B(cartd,nstre,nnode,nunkn,nelem,etype,ptype,ftype);
    dvolu = weigp(igaus)*djacb;
    
    %[chi] = cal_caracteristic_function(igaus,phifunct,element,dim,problembsc,vol_void);
    %chi = 1-vol_void;
    %chi(abs(chi)<1e-12) =  element.material.opt_epsi;

    [Ce] = compute_consitutive_law(element,problembsc,igaus,dim,gamma);
    
    E=zeros(nstre,nnode*nunkn,nelem);                          
    
    for i=1:nstre
        E(i,:,:) = sum(repmat(permute(Ce(i,:,:),[2,1,3]),1,nnode*nunkn,1) .* Bmat,1);
    end
    
    %             inodes=[1 1 2 2 3 3 4 4 ]; icomps=[1 2 1 2 1 2 1 2 ];
    inodes=reshape(repmat([1:nnode],nunkn,1),1,[]); icomps=repmat(1:nunkn,1,nnode);
    %             R=sparse(nndof,nndof);
    
    I_index = zeros( nnode * nunkn * nnode * nunkn * nelem,1);
    J_index = zeros( nnode * nunkn * nnode * nunkn * nelem,1);
    Kij = zeros( nnode * nunkn * nnode * nunkn *  nelem,1);
    k=1;
    % under-diagonal entries & transpose
    for i=1:nnode*nunkn
        ik=nunkn*(inodes(i)-1)+icomps(i);
        it=nunkn*(element.conectivities(:,inodes(i))-1)+icomps(i);
        for j=1:i-1
            jl=nunkn*(inodes(j)-1)+icomps(j);
            jt=nunkn*(element.conectivities(:,inodes(j))-1)+icomps(j);
            k_ij=squeeze(sum(Bmat(:,ik,:) .* E(:,jl,:),1));
            
            I_index(k:k+2*length(it)-1,1) = [ it ; jt];
            J_index(k:k+2*length(it)-1,1) = [ jt ; it];
            Kij(k:k+2*length(it)-1,1) = [ dvolu.*k_ij ; dvolu.*k_ij ];
            k = k + 2*length(it);
        end
    end
    
    % diagonal entries
    for i=1:nnode*nunkn
        ik=nunkn*(inodes(i)-1)+icomps(i);
        it=nunkn*(element.conectivities(:,inodes(i))-1)+icomps(i);
        
        k_ij=squeeze(sum(Bmat(:,ik,:) .* E(:,ik,:),1));
        
        I_index(k:k+length(it)-1,1) =  it ;
        J_index(k:k+length(it)-1,1) =  it ;
        Kij(k:k+length(it)-1,1) =  dvolu.*k_ij ;
        k = k + length(it) ;
    end
    StifMat = StifMat + sparse(I_index,J_index,Kij,nndof,nndof);
    post.chi(igaus,:) = gamma;
end
StifMat = 1/2 * (StifMat + StifMat');
=======
function [StifMat,post] = stiff_OPT_symmetric(dim,element,problembsc,coordinatesn,coordinatesa,gamma)
nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

% dirichlet_data = zeros(nnode,nelem);
idx   = zeros(nnode*nunkn,nelem);
ptype = problembsc.problemtype;
ftype = problembsc.phisical_type;
StifMat = sparse(nndof,nndof);
% StifMat_gauss = repmat({sparse(nndof,nndof)},1,element.ngaus);

dirichlet_data = element.conectivities';
% for i=1:nnode
%     dirichlet_data(i,:)= element.conectivities(:,i);
% end

etype = element.type;
for inode=1:nnode
    for idime=1:nunkn
        idx(nunkn*inode-nunkn+idime,:) = nunkn.*dirichlet_data(inode,:)-nunkn+idime;
    end
end
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
post.chi = zeros(ngaus,nelem);

% tic
for igaus=1:ngaus
     [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [Bmat] = cal_B(cartd,nstre,nnode,nunkn,nelem,etype,ptype,ftype);
    dvolu = weigp(igaus)*djacb;
    
    %[chi] = cal_caracteristic_function(igaus,phifunct,element,dim,problembsc,vol_void);
    %chi = 1-vol_void;
    %chi(abs(chi)<1e-12) =  element.material.opt_epsi;

    [Ce] = compute_consitutive_law(element,problembsc,igaus,dim,gamma);
    
    E=zeros(nstre,nnode*nunkn,nelem);                          
    
    for i=1:nstre
        E(i,:,:) = sum(repmat(permute(Ce(i,:,:),[2,1,3]),1,nnode*nunkn,1) .* Bmat,1);
    end
    
    %             inodes=[1 1 2 2 3 3 4 4 ]; icomps=[1 2 1 2 1 2 1 2 ];
    inodes=reshape(repmat([1:nnode],nunkn,1),1,[]); icomps=repmat(1:nunkn,1,nnode);
    %             R=sparse(nndof,nndof);
    
    I_index = zeros( nnode * nunkn * nnode * nunkn * nelem,1);
    J_index = zeros( nnode * nunkn * nnode * nunkn * nelem,1);
    Kij = zeros( nnode * nunkn * nnode * nunkn *  nelem,1);
    k=1;
    % under-diagonal entries & transpose
    for i=1:nnode*nunkn
        ik=nunkn*(inodes(i)-1)+icomps(i);
        it=nunkn*(element.conectivities(:,inodes(i))-1)+icomps(i);
        for j=1:i-1
            jl=nunkn*(inodes(j)-1)+icomps(j);
            jt=nunkn*(element.conectivities(:,inodes(j))-1)+icomps(j);
            k_ij=squeeze(sum(Bmat(:,ik,:) .* E(:,jl,:),1));
            
            I_index(k:k+2*length(it)-1,1) = [ it ; jt];
            J_index(k:k+2*length(it)-1,1) = [ jt ; it];
            Kij(k:k+2*length(it)-1,1) = [ dvolu.*k_ij ; dvolu.*k_ij ];
            k = k + 2*length(it);
        end
    end
    
    % diagonal entries
    for i=1:nnode*nunkn
        ik=nunkn*(inodes(i)-1)+icomps(i);
        it=nunkn*(element.conectivities(:,inodes(i))-1)+icomps(i);
        
        k_ij=squeeze(sum(Bmat(:,ik,:) .* E(:,ik,:),1));
        
        I_index(k:k+length(it)-1,1) =  it ;
        J_index(k:k+length(it)-1,1) =  it ;
        Kij(k:k+length(it)-1,1) =  dvolu.*k_ij ;
        k = k + length(it) ;
    end
    StifMat = StifMat + sparse(I_index,J_index,Kij,nndof,nndof);
    post.chi(igaus,:) = gamma;
end
StifMat = 1/2 * (StifMat + StifMat');
>>>>>>> refs/remotes/origin/master
end