function [P_operator] = compute_P_operator(Msmooth,element,dim)

ndime=dim.ndime; nelem=dim.nelem; nnode=dim.nnode; 

dirichlet_data=zeros(nnode,nelem);
for inode=1:nnode
    dirichlet_data(inode,:)=element.conectivities(:,inode);
end

T_nodal_2_gauss = sparse(nelem,dim.npnod);
% T_nodal_2_gauss = sparse(zeros(nelem,dim.npnod));
    for inode=1:nnode
       T_nodal_2_gauss = T_nodal_2_gauss + sparse([1:nelem],[dirichlet_data(inode,:)],ones(nelem,1),nelem,dim.npnod);
    end

m = T_nodal_2_gauss*sum(Msmooth,2);
P_operator = diag(m)\T_nodal_2_gauss;       
end