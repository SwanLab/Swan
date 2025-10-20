function dom = TraceEye(A)
    s.operation = @(xV) evaluate(A,xV);
    s.ndimf     = A.ndimf;
    s.mesh      = A.mesh;
    dom         = DomainFunction(s);
end

function trcMat = evaluate(A, xV)
    ndim = A.ndimf(1);
    trc = trace(A).evaluate(xV);
    trcMat = pagemtimes(eye(ndim),trc);
end