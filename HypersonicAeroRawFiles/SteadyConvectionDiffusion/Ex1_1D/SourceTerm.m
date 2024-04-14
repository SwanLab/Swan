function res = SourceTerm(x)
% res = SourceTerm(x)
% Source term for the convection-diffusion equation

global problem

if problem == 1
    res = 0;
elseif problem == 2
    res = 1;
elseif problem == 3
    res = sin(pi*x);
elseif problem == 4
    res = 20*exp(-5*(x-1/8))-10*exp(-5*(x-1/4));
elseif problem == 5
    res = 10*exp(-5*x)-4*exp(-x);
end 