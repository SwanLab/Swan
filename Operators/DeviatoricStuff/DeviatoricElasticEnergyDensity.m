function dom = DeviatoricElasticEnergyDensity(u,mu)
    s.operation = @(xV) evaluate(u,mu,xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u,mu,xV)
    N      = u.ndimf;
    muEval = mu.evaluate(xV);
    nEval  = size(muEval,3);
    e      = AntiVoigt(SymGrad(u));
    D      = Voigt(Deviatoric(e));
    A      = VoigtDeviatorNormMaterial(N,nEval);
    dsE    = DDP(D,DDP(A,D));
    fVR    = 0.5*muEval.*dsE.evaluate(xV);
end