function dom = Shape(u)
    s.operation = @(xV) evaluate(u, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    fVR = u.computeShapeFunctions(xV);
end