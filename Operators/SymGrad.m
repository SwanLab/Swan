function dom = SymGrad(u)
    gradU = Grad(u);
    s.operation = @(xV) evaluate(gradU, xV);
    s.ndimf     = gradU.ndimf;
    dom = DomainFunction(s);
    dom = Voigt(dom);    
end

function symGrad = evaluate(gradU, xV)
    grad  = gradU.evaluate(xV);
    gradT = pagetranspose(grad);
    symGrad = 0.5*(grad + gradT);
end