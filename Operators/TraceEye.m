function dom = TraceEye(A)
    s.operation = @(xV) evaluate(A,xV);
    s.ndimf     = A.ndimf;
    dom         = DomainFunction(s);
end

function trcMat = evaluate(A, xV)
    ndim = sqrt(A.ndimf);
    trc = trace(A).evaluate(xV);
    trcMat = pagemtimes(eye(ndim),trc);
end