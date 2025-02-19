function dom = SymGrad(u)
    gradU = Grad(u);
    s.operation = @(xV) evaluate(gradU, xV);
    s.ndimf     = gradU.ndimf;
    s.mesh      = u.mesh;
    dom = DomainFunction(s);
    dom = Voigt(dom);    
end

function symGrad = evaluate(gradU, xV)
    symGradFun = 0.5*(gradU + gradU');
    symGrad = symGradFun.evaluate(xV);
end