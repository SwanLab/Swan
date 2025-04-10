function dom = ShapeDer(u)
    s.operation = @(xV) evaluate(u, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    dNdx = u.evaluateCartesianDerivatives(xV);
    fVR = dNdx;
    % fVR = reshape(grad, [nDimG*nDimf,nPoints, nElem]);
end