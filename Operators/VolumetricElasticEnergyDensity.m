function dom = VolumetricElasticEnergyDensity(u,kappa)
    s.operation = @(xV) evaluate(u,kappa, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u,kappa, xV)
    divu = Divergence(u);
    dbS  = DDP(kappa,DDP(divu,divu));
    fVR  = 0.5*dbS.evaluate(xV);
end