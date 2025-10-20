function dom = SymGrad(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf     = [u.mesh.ndim,u.ndimf];
    s.mesh      = u.mesh;
    dom = DomainFunction(s);  
end

function symGrad = evaluate(u, xV)
    gradU = Grad(u);
    symGradFun = 0.5*(gradU + gradU');
    symGrad = symGradFun.evaluate(xV);
end