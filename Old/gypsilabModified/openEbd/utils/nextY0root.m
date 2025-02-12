function [ R ] = nextY0root(k)
% First root of Y_0 (bessel function of second kind 0 order) bigger than k
% i.e. R = nextY0root(k) is an approximation of the real z satisfying
% - z >= k
% - Y_0(z) = 0
% - $\forall z', abs(Y_0(z')) < 1e-12 \implies z' >= z$

% Find 3 roots around k. 
rho = besselYroots(k,3); 
% Keep roots larger than k
rho = rho(rho>k);
% Keep the smallest one. 
R = rho(1);

end

