function dom = Grad(u)
    s.operation = @(xV) evaluate(u, xV);
    dom = DomainFunction(s);
end

function df = evaluate(f, xV)
    dNdx = f.evaluateCartesianDerivatives(xV);
    nDofs   = f.nDofs;
    nDimf   = f.ndimf;
    nDimG   = size(dNdx, 1);
    nNodeE  = size(dNdx, 2);
    nPoints = size(dNdx, 3);
    nElem   = size(dNdx, 4);
    connec = f.getConnec();
    fV = reshape(f.fValues', [1 nDofs]);
    df = zeros(nDimG,nDimf, nPoints, nElem);
    for iDimG = 1:nDimG
        for jDimf = 1:nDimf
            for kNodeE = 1:nNodeE
                dNdxIK = squeezeParticular(dNdx(iDimG, kNodeE,:,:),[1 2]);
                iDofE = nDimf*(kNodeE-1)+jDimf;
                dofs = connec(:,iDofE);
                fkj = repmat(fV(dofs),[nPoints 1]);
                dFij = dNdxIK.*fkj;
                df(iDimG,jDimf,:,:) = squeezeParticular(df(iDimG,jDimf,:,:),[1 2]) + dFij;
            end
        end
    end
    df = reshape(df, [nDimG*nDimf,nPoints, nElem]);
end