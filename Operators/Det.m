function dom = Det(A)
    s.operation = @(xV) evaluate(A,xV);
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,xV)
    F = A.evaluate(xV);
    nGaus = size(xV,2);
    nElem = size(F,4);
    jac(1,1,:,:) = MatrixVectorizedInverter.computeDeterminant(F);
    fVR   = reshape(jac, [1 1 nGaus nElem]);
end