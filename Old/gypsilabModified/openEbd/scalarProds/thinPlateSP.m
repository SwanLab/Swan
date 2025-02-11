function[SP] = thinPlateSP(a,rho,R)
% Scalar product with x -> x^2log(Rx)

if length(a)==1
    b= 1;
else
    b = a(2);
    a = a(1);
end

    function[dDg] = dDeltaG(r)
        dDg = 4./r;
    end
    function[DG] = DeltaG(r)
        DG = 4*log(R*r) + 4;
    end
    function[diffG] = dG(r)
        diffG = 2*r.*log(R*r) + r;
    end
    function[res] = ep(rho,r)
        assert(isscalar(r));
        rho = rho(:);
        res = Cp(rho).*besselj(0,rho*r);
    end
    function[res] = dep(rho,r)
        assert(isscalar(r));
        rho = rho(:);
        res = -rho.*Cp(rho).*besselj(1,rho*r);
    end

    function[res] = F(r)
        if r == 0
            res = 2*pi*(-ep(rho,0).*4./rho.^2);
        else
            res = 2*pi*r*(-ep(rho,r).*dDeltaG(r)./rho.^2 ...
                + dep(rho,r).*DeltaG(r)./rho.^2 ...
                + ep(rho,r)*dG(r));
        end
        
    end

SP = F(b) - F(a);

end


