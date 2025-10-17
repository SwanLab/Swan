function dom = Grad(u)
    s.operation = @(xV) u.computeGrad(xV);
    s.ndimf     = [u.mesh.ndim,u.ndimf];
    s.mesh      = u.mesh;
    dom        = DomainFunction(s);
end