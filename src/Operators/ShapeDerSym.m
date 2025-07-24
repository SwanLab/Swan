function dom = ShapeDerSym(u,iDof)
    s.operation = @(xV) evaluate(u,iDof,xV);
    s.mesh      = u.mesh;
    dom         = DomainFunction(s);
end

function fVR = evaluate(u,iDof,xV)
    dN     = ShapeDer(u);
    dNdx   = dN.evaluate(xV);
    nnodeE = u.mesh.nnodeElem;
    ndim   = u.mesh.ndim;
    ndofE  = nnodeE*ndim;
    nGauss = size(xV,2);
    nElem  = u.mesh.nelem;
    gradN  = zeros(ndim,ndim,ndofE,nGauss,nElem);
    for i=1:ndim
        gradN(i,:,i:ndim:end,:,:) = dNdx;
    end
    gradN = squeezeParticular(gradN(:,:,iDof,:,:),3);
    fVR = 0.5*(gradN + pagetranspose(gradN));
end