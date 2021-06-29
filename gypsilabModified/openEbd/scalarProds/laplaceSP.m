function [ out ] = laplaceSP(a,rho)

if length(a)==1
    b= 1;
else
    b = a(2);
    a = a(1);
end
F = @(r)(Cp(rho).*besselj(0,rho*r));
out = 2*pi*(F(b) - F(a));


end

