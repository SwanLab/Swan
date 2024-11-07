function dom = Grad(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf = u.mesh.ndim*u.ndimf;
    dom = DomainFunction(s);
end

function grad = evaluate(u, xV)
    dNdx = u.evaluateCartesianDerivatives(xV);
    uV   = u.getValuesByElem();
    uV   = permute(uV,[1 2 4 3]);    
    grad = pagemtimes(dNdx,uV);
end


