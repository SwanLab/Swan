function dom = ElasticEnergyDensity(C,u)
    s.operation = @(xV) evaluate(C,u,xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(C,u,xV)
    e   = SymGrad(u);
    diE = DDP(e,DDP(C,e));
    fVR = 0.5*diE.evaluate(xV);
end