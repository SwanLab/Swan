function [fg,A_nodal_2_gauss,B_nodal_2_gauss,dvolu] = interpol(fn,element,dim,problembsc,coordinates)
%
ndime=dim.ndime; nelem=dim.nelem; nnode=dim.nnode; 

lnods=zeros(nnode,nelem);
for inode=1:nnode
    lnods(inode,:)=element.conectivities(:,inode);
end
fe=zeros(nnode,nelem);
for inode=1:nnode
    fe(inode,:)=fn(lnods(inode,:));
end
neres=0;
ptype=problembsc.problemtype;
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
fg=zeros(ngaus,nelem);
A_nodal_2_gauss = sparse(nelem,dim.npnod);
B_nodal_2_gauss = sparse(nelem,dim.npnod);
% A_nodal_2_gauss = sparse(zeros(nelem,dim.npnod));
% B_nodal_2_gauss = sparse(zeros(nelem,dim.npnod));
for igaus=1:ngaus
    [~,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinates,coordinates,ptype);
    dvolu = weigp(igaus)*djacb;
    [shape,~,~] = shape_deriv_functions(igaus,posgp,ptype,element.type,nnode,neres);
    for inode=1:nnode
       fg(igaus,:) = fg(igaus,:) + shape(inode)*fe(inode,:);
      
       A_nodal_2_gauss = A_nodal_2_gauss + sparse([1:nelem],[lnods(inode,:)],ones(nelem,1)*shape(inode),nelem,dim.npnod);
      % B_nodal_2_gauss = B_nodal_2_gauss + sparse([1:nelem],[lnods(inode,:)],dvolu*shape(inode),nelem,dim.npnod); 
    end
end
% B_nodal_2_gauss = diag(dvolu)*A_nodal_2_gauss;       
end

