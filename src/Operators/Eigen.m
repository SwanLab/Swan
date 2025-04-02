function dom = Eigen(A)
    s.operation = @(xV) evaluate(A,xV);
    s.mesh = A.mesh;
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,xV)
    F = A.evaluate(xV);
    fVR = pageeig(F);
end