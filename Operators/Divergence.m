function dom = Divergence(u)
    s.operation = @(xV) evaluate(u, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    dNdx = u.evaluateCartesianDerivatives(xV);
    nDofs   = u.nDofs;
    nDimf   = u.ndimf;
    nDimG   = size(dNdx, 1);
    nNodeE  = size(dNdx, 2);
    nPoints = size(dNdx, 3);
    nElem   = size(dNdx, 4);
    connec = u.getDofConnec();
    fV = reshape(u.fValues', [1 nDofs]);
    div = zeros(1, nPoints, nElem);
    for iDim = 1:nDimG
        for kNodeE = 1:nNodeE
            dNdxIK = squeezeParticular(dNdx(iDim, kNodeE,:,:),[1 2]);
            iDofE = nDimf*(kNodeE-1)+iDim;
            dofs = connec(:,iDofE);
            fKJ = repmat(fV(dofs),[nPoints 1]);
            gradIJ= dNdxIK.*fKJ;
            div(1,:,:) = squeezeParticular(div(1,:,:),1) + gradIJ;
        end
    end

    fVR = div;
end