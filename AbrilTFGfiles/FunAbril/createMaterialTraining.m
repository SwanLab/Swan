function material = createMaterialTraining(mesh,nSubdomains,inclusionType,geomFun)
    [young,poisson] = computeElasticProperties(mesh,nSubdomains,inclusionType,geomFun);
    s.type          = 'ISOTROPIC';
    s.ptype         = 'ELASTIC';
    s.ndim          = mesh.ndim;
    s.young         = young;
    s.poisson       = poisson;
    tensor          = Material.create(s);
    material        = tensor;
end

function [young,poisson] = computeElasticProperties(mesh,nSubdomains,inclusionType,geomFun)
    E  = 1;
    nu  = 1/3;
    mD = createMeshDomain(mesh,nSubdomains);

    switch inclusionType
        case {'Hole','HoleRaul'}
            young   = ConstantFunction.create(E,mD);
            poisson = ConstantFunction.create(nu,mD);
        case 'Material'
            E2 = E/1000;

            xmin = min(mesh.coord(:,1));
            xmax = max(mesh.coord(:,1));
            ymin = min(mesh.coord(:,2));
            ymax = max(mesh.coord(:,2));
            Lx = xmax - xmin; % tamaño real en X
            Ly = ymax - ymin; % tamaño real en Y

            fBase             = geomFun;
            fPeriodic = @(x) fBase( cat(1, xmin + mod(x(1,:,:) - xmin, Lx), ...
                                        ymin + mod(x(2,:,:) - ymin, Ly) ));
            
            f = @(x) (fPeriodic(x) <= 0)*E2 + ...
                     (fPeriodic(x) > 0)*E;
    
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