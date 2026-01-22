function material = createMaterialTraining(mesh, r,nSubdomains,inclusionType)
    [young,poisson] = computeElasticProperties(mesh,r,nSubdomains,inclusionType);
    s.type          = 'ISOTROPIC';
    s.ptype         = 'ELASTIC';
    s.ndim          = mesh.ndim;
    s.young         = young;
    s.poisson       = poisson;
    tensor          = Material.create(s);
    material        = tensor;
end

function [young,poisson] = computeElasticProperties(mesh,r,nSubdomains,inclusionType)
    E  = 1;
    nu  = 1/3;
    mD = createMeshDomain(mesh,nSubdomains);

    switch inclusionType
        case {'Hole','HoleRaul'}
            young   = ConstantFunction.create(E,mD);
            poisson = ConstantFunction.create(nu,mD);
        case 'Material'
            E2 = E/1000;
            xmax = max(mesh.coord(:,1));
            ymax = max(mesh.coord(:,2));
            xmin = min(mesh.coord(:,1));
            ymin = min(mesh.coord(:,2));
            Lx = xmax-xmin;
            Ly = ymax-ymin;
    
            f = @(x) ...
                ( sqrt( (mod(x(1,:,:) - xmin, Lx) - Lx/2).^2 + ...
                        (mod(x(2,:,:) - ymin, Ly) - Ly/2).^2 ) < r ) * E2 + ...
                ( sqrt((mod(x(1,:,:) - xmin, Lx) - Lx/2).^2 + ...
                       (mod(x(2,:,:) - ymin, Ly) - Ly/2).^2 ) >= r ) * E;
    
            young   = AnalyticalFunction.create(f, mD);
            poisson = ConstantFunction.create(nu, mD);
   end
end


function mD = createMeshDomain(mR,nSubdomains)
    if sum(nSubdomains > 1)>= 1
        s.nsubdomains   = nSubdomains; %nx ny
        s.meshReference = mR;
        s.tolSameNode   = 1e-10;
        m = MeshCreatorFromRVE2D(s);
        [mD,~,~,~,~,~,~] = m.create();
    else
        mD = mR;
    end
end