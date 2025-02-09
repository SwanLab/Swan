function[C] = Cp(rho)
% C = Cp(rho)
% C_p(rho) = sqrt{\frac{1}{2\pi\int_{0}^1 r rho^2 J_1(rho r)^2 dr}}
% This constant is chosen so that e_p(r) = C_p(rho) J_0(rho_p r) has a unit
% H_0^1(B) norm. (that is \int_{B}|\nabla e_p|^2 = 1).


C = sqrt(1./(pi*rho.*(besselj(0, rho).^2.*rho+besselj(1, rho).^2.*rho-2*besselj(0, rho).*besselj(1, rho))));
C(rho==0) = 1;


end