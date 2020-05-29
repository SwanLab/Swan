function [out] = coeff2func(alpha,rho,x)
% out = coeff2func(alpha,rho,x) 
% Returns the values at x of the function 
% \sum_{p} alpha_p e_p(r) 
% Where e_p(r) = C_p J_0(rho_p r) (see function Cp.m for a definition of 
% C_p). 

y = x(:)';
C = Cp(rho);

J0vals = besselj(0,rho(:)*y);
out = J0vals'*(C(:).*alpha(:));

reshape(out,size(x,1),size(x,2));

end

