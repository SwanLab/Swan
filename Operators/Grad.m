function dom = Grad(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf = u.mesh.ndim*u.ndimf;
    dom = DomainFunction(s);
end

function grad = evaluate(u, xV)
    dNdx = u.evaluateCartesianDerivatives(xV);
    uF   = u.getValuesByElem();
    grad = pagemtimes(dNdx,uF);
end


