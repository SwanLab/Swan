function dom = ShapeDerSym(u,iDof)
    s.operation = @(xV) evaluate(u,iDof,xV);
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end

function fVR = evaluate(u,iDof,xV)

    fVR = 0.5*(gradN + pagetranspose(gradN));
end