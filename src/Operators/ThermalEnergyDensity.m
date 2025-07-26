function dom = ThermalEnergyDensity(kappa,u)
    s.operation = @(xV) evaluate(kappa,u,xV);
    s.mesh = u.mesh;
    dom = DomainFunction(s);
end

function fVR = evaluate(u, kappa, xV)
    thE = times(kappa,DP(Grad(u),Grad(u))); 
    fVR = thE.evaluate(xV);
end