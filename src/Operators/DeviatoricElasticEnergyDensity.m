function dom = DeviatoricElasticEnergyDensity(u,mu)
    s.operation = @(xV) evaluate(u,mu,xV);
    s.mesh = u.mesh;
    dom = DomainFunction(s);
end

function fVR = evaluate(u,mu,xV)
    muEval = mu.evaluate(xV);
    ed     = DeviatoricStrain(u);
    dsE    = 2*DDP(ed,ed);
    fVR    = 0.5*muEval.*dsE.evaluate(xV);
end