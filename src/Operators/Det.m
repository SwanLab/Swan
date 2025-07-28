function dom = Det(A)
    s.operation = @(xV) evaluate(A,xV);
    s.mesh = A.mesh;
    dom         = DomainFunction(s);
end

function fVR = evaluate(A,xV)
    F = A.evaluate(xV);
    fVR(1,:,:) = MatrixVectorizedInverter.computeDeterminant(F);
end