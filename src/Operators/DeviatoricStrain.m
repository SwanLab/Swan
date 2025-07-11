function dom = DeviatoricStrain(u)
    s.operation = @(xV) evaluate(u,xV);
    s.mesh = u.mesh;
    dom         = DomainFunction(s);
end

function fEval = evaluate(u, xV)
    e     = SymGrad(u);
    ed    = Deviatoric(e);
    fEval = ed.evaluate(xV);
end