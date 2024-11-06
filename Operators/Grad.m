function dom = Grad(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf = u.mesh.ndim*u.ndimf;
    dom = DomainFunction(s);
end

function grad = evaluate(u, xV)
    dNdx    = u.evaluateCartesianDerivatives(xV);
    nPoints = size(xV, 2);
    uF = u.getValuesByElem();
    uF = repmat(uF,[1 1 1 nPoints]);
    uF = permute(uF,[1 2 4 3]);
    grad = pagemtimes(dNdx,uF);
end


