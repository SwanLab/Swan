function [ Msmooth,emat ] = mass_matrix(dim,problembsc,element,coordinatesn,coordinatesa)
% 
% Remark : lumped matrix is computed according to the comments described in
% FEM its basis & fund,(Zienkiewicz,taylos,zhu), pp 568, vol I 
% nproc = 1,2,3; Row sum procedure, Diagonal scaling, Quadrature using nodal points
% In the case of CST element the number of gauss used here is 3 !
% In the case of 4-node Quadrilateral, ng = 9 ! 

% job = 1, lumped matrix (for smoothing)
% job = 2, full mass matrix (for smoothing)
job = problembsc.smoothing_proc;

nnode=dim.nnode; nelem=dim.nelem; npnod=dim.npnod; ndime=dim.ndime; nnode=dim.nnode;
dirichlet_data = zeros(nnode,nelem);
ptype = problembsc.problemtype;
etype = element.type;
for i=1:nnode
    dirichlet_data(i,:)= element.conectivities(:,i);
end
switch etype
    case {'TRIANGLE'} % SIMPLICIAL, TRIANG, 3 NODES
        switch nnode
            case 3
                ngaus = 3;
        end
        [posgp,weigp] = trian_gauss_const(ndime,ngaus);
    case 'QUAD' % quadrilateral
        switch nnode
            case 4
                ngaus = 9;
            case {8,9}
                ngaus = 16;
        end
        [posgp,weigp] = quad_gauss_const(ndime,ngaus);
end

emat = zeros(nnode,nnode,nelem);
% compute all the elemental mass matrices
for igaus=1:ngaus
    neres = 0;
    [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [shape,deriv,heslo] = shape_deriv_functions(igaus,posgp,ptype,etype,nnode,neres);
    for inode=1:nnode
        for jnode=1:nnode
            emat(inode,jnode,:)=squeeze(emat(inode,jnode,:)) + weigp(igaus)*shape(inode)*shape(jnode)*djacb(:);
        end
    end
end

if (job==1)
    % lumped mass matrix
    elumped = zeros(nnode,nelem);
    Msmooth = zeros(npnod,1);
    [nproc,coeff] = nprocedure(etype,nnode);
    if (nproc==1)
        for inode=1:nnode
            for jnode=1:nnode
                elumped(inode,:)=elumped(inode,:)+squeeze(emat(inode,jnode,:))';
            end
        end
    elseif (nproc==2)
        for inode=1:nnode
            for jnode=1:nnode
                elumped(inode,:)=elumped(inode,:)+squeeze(emat(inode,jnode,:))';
            end
            elumped(inode,:)=elumped(inode,:)*coeff(inode);
        end        
    end
    for inode=1:nnode
        Msmooth = Msmooth + sparse(dirichlet_data(inode,:),1,elumped(inode,:),npnod,1);
    end
elseif (job==2)
    % full mass matrix
% tstart = tic;
    Msmooth = sparse(npnod,npnod);
        for k=1:nnode
            for l=1:nnode
                vmass = squeeze(emat(k,l,:));
                Msmooth = Msmooth + sparse(dirichlet_data(k,:),dirichlet_data(l,:),vmass,npnod,npnod);
            end
        end
% t1 = toc(tstart)
% tstart = tic;
%     for k=1:nnode
%         for l=1:nnode
%             for e=1:nelem
%                 
%                 Msmooth2(dirichlet_data(k,e),dirichlet_data(l,e))= Msmooth(dirichlet_data(k,e),dirichlet_data(l,e))+emat(k,l,e);
%             end
%         end
%     end
%     t2 = toc(tstart)
% 
% t2/t1
end


end

