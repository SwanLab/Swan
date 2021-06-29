function [ rho,alpha ] = mergeQuads(rho1,alpha1,rho2,alpha2)
%[rho,alpha ] = mergeQuads(rho1,alpha1,rho2,alpha2)
% Computes the vectors alpha and rho such that 
% \sum_{p = 1}^P1 \alpha1_p Cp(rho1)J_0( \rho1_p r) + 
% \sum_{p = 1}^P2 \alpha2_p Cp(rho2)J_0( \rho2_p r) =
% \sum_{p = 1}^P \alpha_p Cp(rho)J_0( \rho_p r)
% See Cp.m for definition of Cp, J_0 is the Bessel function of first kind
% order 0. 
% inputs : 
% - rho1, rho2 : arrays of positive reals (the frequencies). 
% - alpha1, alpha2 : arrays (of the same sizes as rho1 and
% rho2 resepctively) of complex numbers. 


rho1_2 = rho1(~ismember(rho1,rho2));
alpha1_2 = alpha1(~ismember(rho1,rho2));
rho2_1 = rho2(~ismember(rho2,rho1));
alpha2_1 = alpha2(~ismember(rho2,rho1));
[rho12,I1,I2] = intersect(rho1,rho2);
alpha12 = alpha1(I1) + alpha2(I2);

rho_merge = [rho1_2; rho2_1; rho12];
alpha_merge = [alpha1_2;alpha2_1;alpha12];

[rho,I] = sort(rho_merge);
alpha = alpha_merge(I);

end

