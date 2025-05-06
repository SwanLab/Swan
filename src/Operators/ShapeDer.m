function dom = ShapeDer(u)
    s.operation = @(xV) evaluate(u, xV);
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    dNdx = u.evaluateCartesianDerivatives(xV);
    fVR = dNdx;
end