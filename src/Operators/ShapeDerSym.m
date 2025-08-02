function dom = ShapeDerSym(u)
    s.operation = @(xV) evaluate(u, xV);
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    gradN = ShapeDer(u);
    symGradN = 0.5*(gradN + gradN');
    fVR = symGradN.evaluate(xV);
end