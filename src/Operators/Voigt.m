function dom = Voigt(u)
    s.operation = @(xV) evaluate(u,xV);
    switch u.ndimf
        case 4
            s.ndimf = 3;
        case 9
            s.ndimf = 6;
    end
    s.mesh = u.mesh;
    dom = DomainFunction(s);
end

function voigtA = evaluate(u,xV)
    matA = u.evaluate(xV);
    ndim = size(matA,1);
    switch ndim
        case 1
            voigtA(1,:,:) = matA(1,1,:,:);
        case 2
            voigtA = applyVoigt2D(matA);
        case 3
            voigtA = applyVoigt3D(matA);
    end
end
   
function voigtA = applyVoigt2D(matA)
    nPoints = size(matA,3);
    nElem = size(matA,4);
    voigtA = zeros(3,nPoints,nElem);
    
    voigtA(1,:,:) = matA(1,1,:,:); % xx
    voigtA(2,:,:) = matA(2,2,:,:); % yy
    voigtA(3,:,:) = matA(1,2,:,:) + matA(2,1,:,:); % xy
end

function voigtA = applyVoigt3D(matA)
    nPoints = size(matA,3);
    nElem = size(matA,4);
    voigtA = zeros(6,nPoints,nElem);
    
    voigtA(1,:,:) = matA(1,1,:,:); % xx
    voigtA(2,:,:) = matA(2,2,:,:); % yy
    voigtA(3,:,:) = matA(3,3,:,:); % zz
    voigtA(4,:,:) = matA(1,2,:,:) + matA(2,1,:,:); % xy
    voigtA(5,:,:) = matA(1,3,:,:) + matA(3,1,:,:); % xz
    voigtA(6,:,:) = matA(2,3,:,:) + matA(3,2,:,:); % yz
end

    
    
