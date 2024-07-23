function dom = VolumetricElasticEnergyDensity(u,kappa)
    s.operation = @(xV) evaluate(u,kappa, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u,kappa, xV)
    divu = Divergence(u);
    dbE  = DDP(kappa,DDP(divu,divu));
    fVR  = 0.5*dbE.evaluate(xV);
end