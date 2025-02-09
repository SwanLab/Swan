function[SP] = X2kSP(a,rho,k)

if length(a)==1
    b= 1;
else
    b = a(2);
    a = a(1);
end

alpha = @(k)(-4*k.^2./rho.^2);
beta = @(k)(F(k,b) - F(k,a));

SP = rec2(alpha,beta,0,k);


    function[res] = F(k,r)
        res = 2*pi*r*(...
            (2*k)^2./rho.^2*G(k-1,r).*dep(rho,r)...
            + dG(k,r)*ep(rho,r));
    end
    function[res] = G(k,r)
        res = r^(2*k);
    end
    function[res] = dG(k,r)
        res = 2*k*r^(2*k - 1);
    end
    function[res] = ep(rho,r)
        res = Cp(rho).*besselj(0,rho*r);
    end
    function[res] = dep(rho,r)
        res = -Cp(rho).*rho.*besselj(1,rho*r);
    end



end