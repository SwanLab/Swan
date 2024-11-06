function dom = Grad(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf = u.mesh.ndim*u.ndimf;
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    dNdx    = u.evaluateCartesianDerivatives(xV);
    
    
    nDofs   = u.nDofs;
    nDimf   = u.ndimf;
    nDimG   = size(dNdx, 1);
    nNodeE  = size(dNdx, 2);
    nPoints = size(dNdx, 3);
    nElem   = size(dNdx, 4);
    connec = u.getDofConnec();
    fV = reshape(u.fValues', [1 nDofs]);

    %tic
    grad = zeros(nDimG,nDimf, nPoints, nElem);
    for iDimG = 1:nDimG
        for jDimf = 1:nDimf
            for kNodeE = 1:nNodeE
                dNdxIK = squeezeParticular(dNdx(iDimG, kNodeE,:,:),[1 2]);
                iDofE = nDimf*(kNodeE-1)+jDimf;
                dofs = connec(:,iDofE);
                fKJ = repmat(fV(dofs),[nPoints 1]);
                gradIJ= dNdxIK.*fKJ;
                grad(iDimG,jDimf,:,:) = squeezeParticular(grad(iDimG,jDimf,:,:),[1 2]) + gradIJ;
            end
        end
    end
    fVR = reshape(grad, [nDimG*nDimf,nPoints, nElem]);
    %toc
  

    %tic

    uF = u.getValuesByElem();
    uF = repmat(uF,[1 1 1 nPoints]);
    uF = permute(uF,[1 2 4 3]);
    grad = pagemtimes(dNdx,uF);
    fVR = reshape(grad, [nDimG*nDimf,nPoints, nElem]);
    %toc


    %norm(fVR(:) - fVR2(:))/norm(fVR(:))


end


