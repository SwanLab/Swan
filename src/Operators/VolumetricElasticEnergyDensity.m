function dom = VolumetricElasticEnergyDensity(u,kappa)
    s.operation = @(xV) evaluate(u,kappa, xV);
    s.mesh = u.mesh;
    dom = DomainFunction(s);
end

function fVR = evaluate(u,kappa, xV)
    N   = u.ndimf;
    ev  = VolumetricStrain(u);
    dbE = DDP(ev,ev).*kappa;
    fVR = 0.5*N*dbE.evaluate(xV);
end