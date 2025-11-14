function dom = Det(A)
    s.operation = @(xV) evaluate(A,xV);
    s.mesh      = A.mesh;
    s.ndimf     = 1;
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,xV)
    F = A.evaluate(xV);
    nGaus = size(xV,2);
    nElem = size(F,4);
    jac(1,:,:) = MatrixVectorizedInverter.computeDeterminant(F);
    fVR   = reshape(jac, [1 nGaus nElem]);
end