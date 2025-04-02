function dom = Partial(u,dimG)
    s.operation = @(xV) evaluate(u, dimG, xV);
    s.ndimf     = u.ndimf;
    dom = DomainFunction(s);
end

function fVR = evaluate(u, dimG, xV)
    gradU = Grad(u).evaluate(xV);
    fVR = squeezeParticular(gradU(dimG,:,:,:),1);
end