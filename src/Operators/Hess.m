function dom = Hess(u)
    s.operation = @(xV) u.computeHess(xV);
    s.ndimf     = u.mesh.ndim^2*u.ndimf;
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end
