function dom = VolumetricStrain(u,kappa)
    s.operation = @(xV) evaluate(u,kappa, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u,kappa, xV)
    divu = Divergence(u);
    dbS  = DDP(kappa,DDP(divu,divu));
    fVR  = dbS.evaluate(xV);
end