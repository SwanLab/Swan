function dom = Det(A)
    s.operation = @(xV) evaluate(A,xV);
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,xV)
    F = A.evaluate(xV);
    jac(1,1,:,:) = MatrixVectorizedInverter.computeDeterminant(F);
    fVR   = jac;
end