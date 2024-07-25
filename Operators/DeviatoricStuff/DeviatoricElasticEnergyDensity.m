function dom = DeviatoricElasticEnergyDensity(u,mu)
    s.operation = @(xV) evaluate(u,mu,xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u,mu,xV)
    N      = u.ndimf;
    muEval = mu.evaluate(xV);
    nEval  = size(muEval,3);
    ed     = DeviatoricStrain(u);
    A      = VoigtDeviatorNormMaterial(N,nEval);
    dsE    = DDP(ed,DDP(A,ed));
    fVR    = 0.5*muEval.*dsE.evaluate(xV);
end