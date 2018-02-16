function [mu_func,mu] = mu_interp_simp(alpha_minus,alpha_plus,E_minus,E_plus,nu_minus,nu_plus,p)
alpha = sym('alpha','real');
% mu = (E_plus*(alpha - alpha_minus)^p)/(2*(alpha_plus - alpha_minus)^p*(nu_plus + 1)) - (E_minus*((alpha - alpha_minus)^p/(alpha_plus - alpha_minus)^p - 1))/(2*(nu_minus + 1));
mu = (E_plus*(alpha - alpha_minus)^p)/((2*nu_plus + 2)*(alpha_plus - alpha_minus)^p) - (E_minus*((alpha - alpha_minus)^p/(alpha_plus - alpha_minus)^p - 1))/(2*nu_minus + 2);
mu_func = matlabFunction(mu);
end