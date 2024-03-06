function dom = Partial(u,dimG)
    s.operation = @(xV) evaluate(u, dimG, xV);
    s.ndimf     = u.ndimf;
    dom = DomainFunction(s);
end

function fVR = evaluate(u, dimG, xV)
    dNdx    = u.evaluateCartesianDerivatives(xV);
    nDofs   = u.nDofs;
    nDimf   = u.ndimf;
    nNodeE  = size(dNdx, 2);
    nPoints = size(dNdx, 3);
    nElem   = size(dNdx, 4);
    connec  = u.getConnec();
    fV      = reshape(u.fValues', [1 nDofs]);
    grad    = zeros(nDimf, nPoints, nElem);
    for jDimf = 1:nDimf
        for kNodeE = 1:nNodeE
            dNdxIK = squeezeParticular(dNdx(dimG, kNodeE,:,:),[1 2]);
            iDofE = nDimf*(kNodeE-1)+jDimf;
            dofs = connec(:,iDofE);
            fKJ = repmat(fV(dofs),[nPoints 1]);
            gradIJ= dNdxIK.*fKJ;
            grad(jDimf,:,:) = squeezeParticular(grad(jDimf,:,:),1) + gradIJ;
        end
    end

    fVR = reshape(grad, [nDimf,nPoints, nElem]);
end