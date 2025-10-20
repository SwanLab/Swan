function dom = VolumetricElasticEnergyDensity(u,kappa)
    s.operation = @(xV) evaluate(u,kappa, xV);
    s.ndimf     = 1;
    s.mesh      = u.mesh;
    dom = DomainFunction(s);
end

function fVR = evaluate(u,kappa, xV)
    ndim = u.ndimf(1);
    ev   = VolumetricStrain(u);
    dbE  = DDP(ev,ev).*kappa;
    fVR  = 0.5*ndim*dbE.evaluate(xV);
end