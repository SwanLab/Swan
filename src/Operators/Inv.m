function dom = Inv(A)
    s.operation = @(xV) evaluate(A,xV);
    s.mesh      = A.mesh;
    s.ndimf     = A.ndimf;
    dom         = DomainFunction(s);
end

function invF = evaluate(A,xV)
    F = A.evaluate(xV);
    invF = MatrixVectorizedInverter.computeInverse(F);
end