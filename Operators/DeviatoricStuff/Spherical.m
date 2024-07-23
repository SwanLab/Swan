function dom = Spherical(A)
    s.operation = @(xV) evaluate(A,xV);
    s.ndimf     = A.ndimf;
    dom         = DomainFunction(s);
end

function Asph = evaluate(A, xV)
    ndim = sqrt(A.ndimf);
    Asph = (1/ndim).*TraceEye(A).evaluate(xV);
end