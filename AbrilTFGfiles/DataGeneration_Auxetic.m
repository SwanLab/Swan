%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale'
p.Sampling   = 'Oversampling';  %'Isolated'/'Oversampling'

%% DATA GENERATION
mR              = createReferenceMesh();
switch p.Training
    case 'Multiscale'
        material          = createMaterial(mR,[1 1]);
        mesh              = mR;
        s.type             = 'discontinuous';
        switch s.type
            case 'continuous'
                bMesh        = mesh.createSingleBoundaryMesh();
                s.lambdaFun  = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
            case 'discontinuous'
                bMesh        = mesh.createBoundaryMesh();
                for i=1:numel(bMesh)
                    s.lambdaFun{i}=LagrangianFunction.create(bMesh{i}.mesh,mesh.ndim, 'P1');
                end
        end
        s.mesh            = mesh;
        s.bMesh           = bMesh;
        s.uFun            = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
        s.material        = material;
        s.dirichletFun    = createDirichletFunction(bMesh,s.type);
        e                 = ElasticHarmonicExtension(s);
        [T,lambda,K,Kcoarse] = e.solve();
   
    case {'EIFEM','EIFisol'}
        [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
        [material]     = createMaterial(mR,nS);
        s.mesh           = mR;
        s.material       = material;
        s.domainIndices  = dI;
        s.nSubdomains    = nS;
        s.type           = 'discontinuous';
        m= EIFEMTraining(s);
        data          = m.train();
        [data.material] = createMaterial(mR,[1 1]);
        data.dirac=true;
        z = OfflineDataProcessor(data);

        EIFEoper = z.computeROMbasis();
        T        = EIFEoper.U;
        mesh     = data.mesh;
        Kcoarse  = EIFEoper.Kcoarse;
end

string = "CellEIFEMtraining.mat";

% Guarda el .mat per cert radi
FileName=fullfile('AbrilTFGfiles','Data',"Auxetic",string);

switch p.Training
    case 'Multiscale'
        save(FileName,"T","Kcoarse");
    case 'EIFEM'
        save(FileName, "EIFEoper","T","Kcoarse","mesh");
end



%% FUNCTIONS

function mS = createReferenceMesh()

    filename = 'DEF_Q4auxL_1.mat';
    % filename = 'DEF_Q4porL_1.mat';
    load(filename,'EIFEoper');
    s.coord    = EIFEoper.MESH.COOR;
    s.connec   = EIFEoper.MESH.CN;

    obj.xmin = min(s.coord(:,1));
    obj.xmax = max(s.coord(:,1));
    obj.ymin = min(s.coord(:,2));
    obj.ymax = max(s.coord(:,2));

    delta = 1e-9;
    s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
    s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[delta,-delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[delta,delta];
    mS = Mesh.create(s);

end

function dF = createDirichletFunction(bMesh,type)
s.mesh = bMesh;
s.order= 1;
s.type = type;
cf = CoarseFunctions(s);
dF = cf.getAnalytical();
end	


function [nS,dI] = defineNumberOfSubdomains(type)
    switch type
        case 'Isolated'
            nS = [1 1]; %nx ny
            dI = [1 1];
        case 'Oversampling'
            nS = [5 5]; %nx ny
            dI = [3 3];
    end
end


function [material] = createMaterial(mesh,nS)
    mD           = createMeshDomain(mesh,nS);
    E       = 1;
    nu      = 1/3;
    young   = ConstantFunction.create(E,mD);
    poisson = ConstantFunction.create(nu,mD);
    s.type          = 'ISOTROPIC';
    s.ptype         = 'ELASTIC';
    s.ndim          = mesh.ndim;
    s.young         = young;
    s.poisson       = poisson;
    material        = Material.create(s);
end


function meshDom=createMeshDomain(mR,nS)
    if sum(nS > 1)>= 1
        s.nsubdomains   = nS; %nx ny
        s.meshReference = mR;
        s.tolSameNode   = 1e-11;
        m = MeshCreatorFromRVE.create(s);
        [meshDom,~,~,~,~,~,~] = m.create();
    else
        meshDom = mR;
    end
end