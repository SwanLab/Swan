function [f_g] = scl_product(nelem,nnode,lnods,emass,f,g)    
% compute intgral_omega(f*g)
% f,g are suposed to be nodal functions

prod = zeros(1,nelem);
ef = zeros(nnode,nelem);
eg = zeros(nnode,nelem);

for i=1:nnode
    ef(i,:)= f(lnods(i,:));
    eg(i,:)= g(lnods(i,:));
end

for i=1:nnode
    for j=1:nnode
        vmass = squeeze(emass(i,j,:))';
        prod = prod + eg(i,:).*vmass.*ef(j,:);
    end
end
f_g = sum(prod);


end

