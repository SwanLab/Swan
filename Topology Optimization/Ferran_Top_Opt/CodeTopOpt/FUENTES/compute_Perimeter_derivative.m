function per_der = compute_Perimeter_derivative(dim,element,epsilon,Stiff,Mass,gamma)

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

% Compute error regularization
p_epsilon = (epsilon^2*Stiff + Mass)\(Mass*(vepsilon - txi));

% Compute per_der
per_der = 4/epsilon*(ones(size(p_epsilon)) + 2*(p_epsilon - vepsilon));

end