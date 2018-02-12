function Perimeter = compute_Perimeter(dim,element,emass,epsilon,Stiff,Mass,vol,gamma)

nelem=dim.nelem;  nnode=dim.nnode;
lnods = zeros(nnode,nelem);
for i=1:nnode
    lnods(i,:)= element.conectivities(:,i);
end

% Computation of caracteristic function
txi = ones(size(gamma));
txi(gamma>0) = 0;

% Compute first regularization
vepsilon = (epsilon^2*Stiff + Mass)\(Mass*txi);

% Compute Per
Per = scl_product(nelem,nnode,lnods,emass,vepsilon,vepsilon-2*txi);

Perimeter = 4/epsilon*(Per + vol);

end