function dom = Deviatoric(A)
    s.operation = @(xV) evaluate(A,xV);
    s.ndimf     = A.ndimf;
    s.mesh      = A.mesh;
    dom         = DomainFunction(s);
end

function Adev = evaluate(A, xV)
    Avol = Spherical(A);
    Adev = A.evaluate(xV) - Avol.evaluate(xV);
end