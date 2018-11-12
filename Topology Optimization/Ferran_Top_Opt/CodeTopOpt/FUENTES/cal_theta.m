function [theta,norm_g,g_ortho,norm_g_ortho,g_phi] = cal_theta(dim,element,g,phi,emass)  
%                       (    (g,phi)   ) 
% compute theta = arccos(--------------)
%                       ( ||g|| ||phi||)
% 
% It is suposed that g and phi are nodal functions. Therefore,
% the appropiated gauss point number should be used. 
%

nelem=dim.nelem;  nnode=dim.nnode;
dirichlet_data = zeros(nnode,nelem);
for i=1:nnode
    dirichlet_data(i,:)= element.conectivities(:,i);
end
[g_phi] = scl_product(nelem,nnode,dirichlet_data,emass,g,phi);   
[g_g] = scl_product(nelem,nnode,dirichlet_data,emass,g,g);   
[phi_phi] = scl_product(nelem,nnode,dirichlet_data,emass,phi,phi);   

norm_g = sqrt(g_g);
norm_phi = sqrt(phi_phi); % verificar que efectivamente es norma uno y borrar
theta = real(acos(g_phi/(norm_g*norm_phi)));

g_ortho = g - g_phi/phi_phi*phi; 
norm_g_ortho = sqrt(scl_product(nelem,nnode,dirichlet_data,emass,g_ortho,g_ortho));
end

