function dom = SymGrad(u)
    s.operation = @(xV) evaluate(u, xV);
    s.ndimf     = obtainNDimFieldsVoigt(u);
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    grad = Grad(u).evaluate(xV);
    nDimf = u.ndimf;
    nDims = size(grad, 1)/nDimf;
    nGaus = size(grad, 2);
    nElem = size(grad, 3);

    gradReshp = reshape(grad, [nDims,nDimf,nGaus,nElem]);
    gradT = permute(gradReshp, [2 1 3 4]);
    symGrad = 0.5*(gradReshp + gradT);
    
    symGrad =  reshape(symGrad, [nDims*nDimf,1, nGaus,nElem]);
    fVR = obtainVoigtFormat(symGrad);
end

function newObj = obtainVoigtFormat(sgr)
    ndim = size(sgr,1);
    switch ndim
        case 4
            newObj = applyVoigt2D(sgr);
        case 9
            newObj = applyVoigt3D(sgr);
    end
end

function ndimfVoigt = obtainNDimFieldsVoigt(u)
    ndimg = u.mesh.ndim;
    ndimf = u.ndimf;
    n     = ndimg*ndimf;
    switch n
        case 4
            ndimfVoigt = 3;
        case 9
            ndimfVoigt = 6;
    end
end

function fV = applyVoigt2D(fV)
    nGaus = size(fV,3);
    nElem = size(fV,4);
    fVal(1,:,:) = fV(1,:,:); % xx
    fVal(2,:,:) = fV(4,:,:); % yy
    fVal(3,:,:) = fV(2,:,:) + fV(3,:,:); % xy
    fV = reshape(fVal, [3 nGaus nElem]);
%             newObj = FGaussDiscontinuousFunction.create(fV,obj.mesh,obj.quadrature);
end

function fV = applyVoigt3D(fV)
    nGaus = size(fV,3);
    nElem = size(fV,4);
    fVal(1,:,:) = fV(1,:,:); % xx
    fVal(2,:,:) = fV(5,:,:); % yy
    fVal(3,:,:) = fV(9,:,:); % zz
    fVal(4,:,:) = fV(2,:,:) + fV(4,:,:); % xy
    fVal(5,:,:) = fV(3,:,:) + fV(7,:,:); % xz
    fVal(6,:,:) = fV(6,:,:) + fV(8,:,:); % yz
    fV = reshape(fVal, [6 nGaus nElem]);
%             newObj = FGaussDiscontinuousFunction.create(fV,obj.mesh,obj.quadrature);
end