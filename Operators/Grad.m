function dom = Grad(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf = u.mesh.ndim*u.ndimf;
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    dNdx    = u.evaluateCartesianDerivatives(xV);
    nDimf   = u.ndimf;
    nDimG   = size(dNdx, 1);
    nPoints = size(dNdx, 3);
    nElem   = size(dNdx, 4);
    uF = u.getValuesByElem();
    uF = repmat(uF,[1 1 1 nPoints]);
    uF = permute(uF,[1 2 4 3]);
    grad = pagemtimes(dNdx,uF);
    fVR = reshape(grad, [nDimG*nDimf,nPoints, nElem]);
end


