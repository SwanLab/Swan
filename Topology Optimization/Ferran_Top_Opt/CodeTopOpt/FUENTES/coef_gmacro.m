function [h1,h2] = coef_gmacro(gamma,mu,lambda)

coeff = (lambda+3*mu)/(lambda+mu);

c1 = (gamma-1)/(coeff*gamma+1)*(coeff+1)/2;
c2 = (gamma-1)*(coeff-2)/(coeff + 2*gamma -1);
h1 = 2*c1;
h2 = c1*c2;


end