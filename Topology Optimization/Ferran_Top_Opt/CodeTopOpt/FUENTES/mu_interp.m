function [mu_func,mu] = mu_interp(alpha_minus,alpha_plus,E_minus,E_plus,nu_minus,nu_plus)
syms c1 c2 c3 c4
syms gamm
% eq = [ (c1*gamma_plus^2 + c2*gamma_plus + 1)/(c4 + c3*gamma_plus)- E_plus/(2*nu_plus + 2);
%     (c1*gamma_minus^2 + c2*gamma_minus + 1)/(c4 + c3*gamma_minus) - E_minus/(2*nu_minus + 2);
%     (c2 + 2*c1*gamma_plus)/(c4 + c3*gamma_plus) - (c3*(c1*gamma_plus^2 + c2*gamma_plus + 1))/(c4 + c3*gamma_plus)^2 + (2*E_plus*(E_minus - E_plus + E_minus*nu_plus - E_plus*nu_minus))/((gamma_plus - gamma_minus)*(nu_plus + 1)^2*(3*E_minus + E_plus - E_minus*nu_plus + E_plus*nu_minus));
%     (c2 + 2*c1*gamma_minus)/(c4 + c3*gamma_minus) - (c3*(c1*gamma_minus^2 + c2*gamma_minus + 1))/(c4 + c3*gamma_minus)^2 + (2*E_minus*(E_minus - E_plus + E_minus*nu_plus - E_plus*nu_minus))/((gamma_plus - gamma_minus)*(nu_minus + 1)^2*(E_minus + 3*E_plus + E_minus*nu_plus - E_plus*nu_minus))];

eq = [ (c1*alpha_plus^2 + c2*alpha_plus + 1)/(c4 + alpha_plus*c3) - E_plus/(2*nu_plus + 2);
       (c1*alpha_minus^2 + c2*alpha_minus + 1)/(c4 + alpha_minus*c3) - E_minus/(2*nu_minus + 2);
       (c2 + 2*alpha_plus*c1)/(c4 + alpha_plus*c3) - (c3*(c1*alpha_plus^2 + c2*alpha_plus + 1))/(c4 + alpha_plus*c3)^2 + (2*E_plus*(E_minus - E_plus + E_minus*nu_plus - E_plus*nu_minus))/((alpha_plus - alpha_minus)*(nu_plus + 1)^2*(3*E_minus + E_plus - E_minus*nu_plus + E_plus*nu_minus));
(c2 + 2*alpha_minus*c1)/(c4 + alpha_minus*c3) - (c3*(c1*alpha_minus^2 + c2*alpha_minus + 1))/(c4 + alpha_minus*c3)^2 + (2*E_minus*(E_minus - E_plus + E_minus*nu_plus - E_plus*nu_minus))/((alpha_plus - alpha_minus)*(nu_minus + 1)^2*(E_minus + 3*E_plus + E_minus*nu_plus - E_plus*nu_minus))];
 

coef = (struct2array(solve(eq,[c1,c2,c3,c4])));
mu = simplify((coef(1)*gamm^2 + coef(2)*gamm + 1)/(coef(3)*gamm + coef(4)));
mu_func = matlabFunction(mu);
end