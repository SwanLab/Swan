function dom = VolumetricStrain(u)
    s.operation = @(xV) evaluate(u,xV);
    dom         = DomainFunction(s);
end

function fEval = evaluate(u, xV)
    e     = AntiVoigt(SymGrad(u));
    ev    = Voigt(Spherical(e));
    fEval = ev.evaluate(xV);
end