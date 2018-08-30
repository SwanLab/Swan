function P = derivative_constitutive_tensor(element,dim,gamma)

c1 = element.polarization_sym_part_fourth_order_tensor;
c2 = element.polarization_doble_product_second_identity_tensor;

if nargin(c1) == 0
   c1_gamma = 0;
else 
    c1_gamma = c1(gamma);
 end

if nargin(c2) == 0
   c2_gamma = 0;
else
    c2_gamma = c2(gamma);
 end


P = zeros(dim.nstre,dim.nstre,length(gamma));

P(1,1,:) = c1_gamma + c2_gamma;
P(2,2,:) = c1_gamma + c2_gamma;
P(1,2,:) = c2_gamma;
P(2,1,:) = c2_gamma;
P(3,3,:) = c1_gamma;


end