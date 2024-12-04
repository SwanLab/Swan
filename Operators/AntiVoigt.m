function dom = AntiVoigt(A)
    s.operation = @(xV) evaluate(A,xV);
    switch A.ndimf
        case 3
            s.ndimf = 4;
        case 6
            s.ndimf = 9;
    end
    s.mesh = A.mesh;
    dom = DomainFunction(s);
end

function AntiVoigt = evaluate(A,xV)
    vecA = A.evaluate(xV);
    ndim = size(vecA,1);
    switch ndim
        case 3
            AntiVoigt = applyAntiVoigt2D(vecA);
        case 6
            AntiVoigt = applyAntiVoigt3D(vecA);
    end
end

function AntiVoigtA = applyAntiVoigt2D(vecA)
    nPoints = size(vecA,2);
    nElem = size(vecA,3);
    AntiVoigtA = zeros(2,2,nPoints,nElem);

    AntiVoigtA(1,1,:,:) = vecA(1,:,:);
    AntiVoigtA(1,2,:,:) = vecA(3,:,:)./2;
    AntiVoigtA(2,1,:,:) = vecA(3,:,:)./2;
    AntiVoigtA(2,2,:,:) = vecA(2,:,:);
end

function AntiVoigtA = applyAntiVoigt3D(vecA)
    nPoints = size(vecA,2);
    nElem = size(vecA,3);
    AntiVoigtA = zeros(3,3,nPoints,nElem);

    AntiVoigtA(1,1,:,:) = vecA(1,:,:);
    AntiVoigtA(1,2,:,:) = vecA(4,:,:)./2;
    AntiVoigtA(1,3,:,:) = vecA(5,:,:)./2;
    AntiVoigtA(2,1,:,:) = vecA(4,:,:)./2;
    AntiVoigtA(2,2,:,:) = vecA(2,:,:);
    AntiVoigtA(2,3,:,:) = vecA(6,:,:)./2;
    AntiVoigtA(3,1,:,:) = vecA(5,:,:)./2;
    AntiVoigtA(3,2,:,:) = vecA(6,:,:)./2;
    AntiVoigtA(3,3,:,:) = vecA(3,:,:);
end