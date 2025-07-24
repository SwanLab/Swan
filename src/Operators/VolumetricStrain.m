function dom = VolumetricStrain(u)
    s.operation = @(xV) evaluate(u,xV);
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end

function fEval = evaluate(u, xV)
    e     = SymGrad(u);
    ev    = Spherical(e);
    fEval = ev.evaluate(xV);
end