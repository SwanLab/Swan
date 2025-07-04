function dom = ShapeDer(u)
    s.operation = @(xV) evaluate(u, xV);
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    fVR = u.evaluateCartesianDerivatives(xV);
end