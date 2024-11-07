function dom = SymGrad(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf = u.ndimf*u.mesh.ndim;
    dom = DomainFunction(s);
    dom = Voigt(dom);    
end

function symGrad = evaluate(u, xV)
    grad = Grad(u).evaluate(xV);
    gradT = pagetranspose(grad);
    symGrad = 0.5*(grad + gradT);
end