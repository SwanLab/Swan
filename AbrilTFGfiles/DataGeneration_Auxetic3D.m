%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

% l=1.6:0.1:1.8;
l=1.8;
p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale', ('EIFisol')
p.Sampling   = 'Oversampling';  %'Isolated'/'Oversampling'
p.Inclusion  = 'Hole';        % 'Material'/'Hole'
p.nelem      = 15;
meshName     = p.nelem+"x"+p.nelem;

%% DATA GENERATION

mR              = createReferenceMesh(p);
switch p.Training
    case 'Multiscale'
        p.Sampling = 'Isolated';
        [material, mT]   = createMaterial(mR,[1 1],'Material',g);
        mesh       = mR;
        s.type             = 'discontinuous';
        switch s.type
            case 'continuous'
                bMesh      = mesh.createSingleBoundaryMesh();
                s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
            case 'discontinuous'
                bMesh        = mesh.createBoundaryMesh();
                for i=1:numel(bMesh)
                    s.lambdaFun{i}=LagrangianFunction.create(bMesh{i}.mesh,mesh.ndim, 'P1');
                end
        end
        s.mesh          = mesh;
        s.bMesh         = bMesh;
        s.material      = material;
        s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
        s.dirichletFun  = createDirichletFunction(bMesh,s.type);
        e  = ElasticHarmonicExtension(s);
        [T,lambda,K,Kcoarse] = e.solve();
        V  = Integrator.compute(mT.designVariable.fun,mR,2);

    case {'EIFEM','EIFisol'}
        [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
        [material,mT]    = createMaterial(mR,nS,'Material',g);
        s.mesh           = mR;
        s.material       = material;
        s.domainIndices  = dI;
        s.nSubdomains    = nS;
        s.type           = 'continuous';
        s.Coarseorder    = 2;
        m= EIFEMTraining(s);
        data              = m.train();
        data.dirac        = true;
        data.type         = s.type;
        data.Coarseorder  = 2;
        data.material     = createMaterial(mR,[1 1],'Material',g);
        z = OfflineDataProcessor(data);

        EIFEoper = z.computeROMbasis();
        T        = EIFEoper.U;
        mesh     = data.mesh;
        Kcoarse  = EIFEoper.Kcoarse;
        vT       = mR.computeVolume();
end


string = "Auxetic3DIsolated";

% Guarda el .mat per cert radi
FileName=fullfile('AbrilTFGfiles','Data',"Cube",p.Training,meshName,string);
% FileName=fullfile('AbrilTFGfiles','Data',"Cube",'EIFisol',meshName,string);

switch p.Training
    case 'Multiscale'
        save(FileName,"T","Kcoarse","mesh","R");
    case {'EIFEM','EIFisol'}
        save(FileName, "EIFEoper","T","Kcoarse","mesh");
end
    


%% FUNCTIONS

function mS = createReferenceMesh(p)
    switch p.Inclusion
        case 'Material'
            mS = createStructuredMesh(p);
        case 'Hole'
            mS = createStructuredMesh(p);
            lvSet    = createLevelSetFunction(mS);
            uMesh    = computeUnfittedMesh(mS,lvSet);
            mR       = uMesh.createInnerMesh();
            mS       = SmoothMesh(mR);
    end
end

function createStructuredMesh(p)
    n =p.nelem;
    m = TetraMesh(1.5,1,1,n,n,n);

    s.coord=m.coord;
    s.connec=m.connec;
    maxC= max(s.coord);
    minC = min(s.coord);
    
    obj.xmin = minC(1);
    obj.xmax = maxC(1);
    obj.ymin = minC(2);
    obj.ymax = maxC(2);
    obj.zmin = minC(3);
    obj.zmax = maxC(3);

    delta=1E-9;
    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[delta,0,0];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)-[delta,0,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)+[delta,0,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[delta,0,0];

    s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

    s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

    s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,delta];

    s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,delta];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];

    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==maxC(2),:)-[0,delta,0];

    s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,2)==minC(2),:)+[0,delta,0];

    mS = Mesh.create(s);
end

function levelSet = createLevelSetFunction(bgMesh)
    gPar.type           = 'Auxetic3D';
    gPar.length         = 2;
    gPar.height         = 2;
    gPar.xCoorCenter    = 0;
    gPar.yCoorCenter    = 0;
    gPar.theta          = 58;
    gPar.beta           = 64;
    gPar.thickness      = 0.3;
    
    g         = GeometricalFunction(gPar);
    phiFun    = g.computeLevelSetFunction(bgMesh);
    levelSet  = phiFun.fValues;
end


function dF = createDirichletFunction(bMesh)
s.mesh = bMesh;
s.type = 'continuous';
cf = CoarseFunctions(s);
dF = cf.getAnalytical();
end	


function [nS,dI] = defineNumberOfSubdomains(type)
    switch type
        case 'Isolated'
            nS = [1 1 1]; %nx ny
            dI = [1 1 1];
        case 'Oversampling'
            nS = [3 3 3]; %nx ny
            dI = [2 2 2];
    end
end


function [material,m] = createMaterial(mesh,nSubdomains,g)
    mD = createMeshDomain(nSubdomains,mesh);
    s.mesh          = mD;
    s.geomFun       = g;
    m = MaterialTraining(s);
    material = m.create();
end

function mD = createMeshDomain(nS,mesh)
    if sum(nS > 1)>= 1
        s.nsubdomains   = nS; %nx ny
        s.meshReference = mesh;
        s.tolSameNode   = 1e-11;
        m = MeshCreatorFromRVE.create(s);
        [mD,~,~,~,~,~,~] = m.create();
    else
        mD = mesh;
    end
end