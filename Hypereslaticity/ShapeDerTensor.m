function dom = ShapeDerTensor(u)
    s.operation = @(xV) evaluate(u, xV);
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end

function gradN = evaluate(u, xV)
    dNdx   = u.evaluateCartesianDerivatives(xV);
    nnodeE = u.mesh.nnodeElem;
    ndim   = u.mesh.ndim;
    ndofE  = nnodeE*ndim;
    nGauss = size(xV,2);
    nElem  = u.mesh.nelem;
    gradN  = zeros(ndim,ndim,ndofE,nGauss,nElem);
    for i=1:ndim
        gradN(i,:,i:ndim:end,:,:) = dNdx;
    end
end