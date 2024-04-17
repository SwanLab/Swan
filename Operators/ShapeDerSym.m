function dom = ShapeDerSym(u)
    s.operation = @(xV) evaluate(u, xV);
    dom = DomainFunction(s);
end

function fVR = evaluate(u, xV)
    dNdx = u.evaluateCartesianDerivatives(xV);
    B = computeB(u,dNdx);
    fVR = B;
end

function B = computeB(u,dNdx)
    ndimf = u.ndimf;
    switch ndimf
        case 1
            B = computeBin1D(u,dNdx);
        case 2
            B = computeBin2D(u,dNdx);
        case 3
            B = computeBin3D(u,dNdx);
    end
end

function [B] = computeBin1D(u,deriv)
    nVoigt = 2;
    nDimf = u.ndimf;
    nNode = size(deriv,2);
    nElem = size(u,4);
    nDofs = nDimf*nNode;
    B = zeros(nVoigt,nDofs,nGaus,nElem);
    for inode = 1:nNode
        j = nDimf*(inode-1) + 1;
        B(1,j,:,:) = deriv(1,inode,:,:);
        B(2,j,:,:) = deriv(2,inode,:,:);
    end
end

function B = computeBin2D(u,deriv)
    nVoigt = 3;
    nDimf = u.ndimf;
    nNodE = size(deriv,2);
    nDofE = nNodE*nDimf;
    nGaus = size(deriv,3);
    nElem = size(deriv,4);
    B = zeros(nVoigt,nDofE,nGaus,nElem);
    for iNode = 1:nNodE
        j = nDimf*(iNode-1)+1;
        B(1,j,:,:)   = deriv(1,iNode,:,:);
        B(2,j+1,:,:) = deriv(2,iNode,:,:);
        B(3,j,:,:)   = deriv(2,iNode,:,:);
        B(3,j+1,:,:) = deriv(1,iNode,:,:);
    end
end

function B = computeBin3D(u,deriv)
    nVoigt = 6;
    nNode = size(deriv,2);
    nGaus = size(deriv,3);
    nElem = size(deriv,4);
    B = zeros(nVoigt,nNode,nGaus,nElem);
    for inode = 1:nNode
        j = u.ndimf*(inode-1)+1;
        % associated to normal strains
        B(1,j,:,:)   = deriv(1,inode,:,:);
        B(2,j+1,:,:) = deriv(2,inode,:,:);
        B(3,j+2,:,:) = deriv(3,inode,:,:);
        % associated to shear strain, gamma12
        B(4,j,:,:)   = deriv(2,inode,:,:);
        B(4,j+1,:,:) = deriv(1,inode,:,:);
        % associated to shear strain, gamma13
        B(5,j,:,:)   = deriv(3,inode,:,:);
        B(5,j+2,:,:) = deriv(1,inode,:,:);
        % associated to shear strain, gamma23
        B(6,j+1,:,:) = deriv(3,inode,:,:);
        B(6,j+2,:,:) = deriv(2,inode,:,:);
    end
end