function [f_g] = scl_product(nelem,nnode,dirichlet_data,emass,f,g)    
% compute intgral_omega(f*g)
% f,g are suposed to be nodal functions

prod = zeros(1,nelem);
ef = zeros(nnode,nelem);
eg = zeros(nnode,nelem);

for i=1:nnode
    ef(i,:)= f(dirichlet_data(i,:));
    eg(i,:)= g(dirichlet_data(i,:));
end

for i=1:nnode
    for j=1:nnode
        vmass = squeeze(emass(i,j,:))';
        prod = prod + eg(i,:).*vmass.*ef(j,:);
    end
end
f_g = sum(prod);


end

