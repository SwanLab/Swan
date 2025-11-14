function dom = Grad(u)
    s.operation = @(xV) u.computeGrad(xV);
    s.ndimf     = [u.ndimf,u.mesh.ndim];
    s.mesh      = u.mesh;
    dom        = DomainFunction(s);
end