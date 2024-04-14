function res = u(x,a,nu)
% res = u(x,a,nu)
% Analytical solution of a one-dimensional convection-diffusion problem
% with essential boundary conditions on both ends.
%  

global problem

if problem == 1
    res = (1-exp(x*a/nu))/(1-exp(a/nu));
elseif problem == 2
    res = (x + (1 - exp(a/nu*x))/((exp(a/nu)-1)))/a;
elseif problem == 3
    aux = pi*(a^2+nu^2*pi^2);
    e = exp(a/nu);
    c1 = (-aux+a*(e+1))/(aux*(e-1));
    c2 = (aux-2*a)/(aux*(e-1));
    res = c1 + c2*exp(a*x/nu) + nu*pi*(sin(pi*x)-a*cos(pi*x)/(nu*pi))/aux;
elseif problem == 4
    res = 2*exp(-5*x+5/8)*(-2+exp(5/8))/(5*nu+a) - ...     
        (5*nu*exp(5/8)^7-2*exp(5/8) - 4*exp(a/nu+5) + 2*exp(a/nu+45/8) + 4 + a*exp(5/8)^7)/...
        ((5*nu+a)*exp(5/8)^7*(exp(a/nu)-1)) + ...
        (4-2*exp(5/8)+5*nu*exp(5/8)^7+a*exp(5/8)^7-4*exp(5/8)^8+2*exp(5/8)^9)*exp(a*x/nu)/...
        (exp(5/8)^7*(5*exp(a/nu)*nu+exp(a/nu)*a-5*nu-a));
elseif problem == 5
    aux = 5*nu^2 + 6*a*nu + a^2;
    res = ( exp(a/nu)*(-18*nu-2*a)+...
        exp(a*x/nu)*(5*nu^2 + 2*a + a^2 + 6*nu*a + 18*nu) - aux + ...
        2*(exp(-5*x) - exp(-5) + exp(-(-a*x+5*nu)/nu) -exp(-(5*x*nu-a)/nu))*(nu + a) - ...    
        4*(exp(-x) + exp(-(-a*x+nu)/nu) - exp(-(x*nu-a)/nu) - exp(-1))*(5*nu + a))/...    
        (aux*(exp(a/nu)-1));
end
