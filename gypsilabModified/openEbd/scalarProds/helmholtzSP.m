function [ out ] = helmholtzSP( a,rho,k )

if length(a)==1
    b = 1;
else
    b = a(2);
    a = a(1);
end

C = Cp(rho);
F = @(r)(C.*k.*r.*(bessely(1, k*r).*besselj(0, rho*r).*rho-...
    besselj(1, rho*r).*bessely(0, k*r)*k).*rho./(k^2-rho.^2));
out = 2*pi*(F(b) - F(a + (a==0)*10*eps));



end

