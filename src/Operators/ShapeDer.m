function dom = ShapeDer(u)
    s.operation = @(xV) evaluate(u, xV);
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end

function gradN = evaluate(u, xV)
    dNdx   = u.evaluateCartesianDerivatives(xV);
    nnodeE = u.mesh.nnodeElem;
    ndim   = u.mesh.ndim;
    ndimf  = u.ndimf;
    ndofE  = nnodeE*ndimf;
    nGauss = size(xV,2);
    nElem  = u.mesh.nelem;
    gradN  = zeros(ndimf,ndim,ndofE,nGauss,nElem);
    for i=1:ndimf
        gradN(i,:,i:ndimf:end,:,:) = dNdx;
    end
end