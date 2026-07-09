%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

% l=1.6:0.1:1.8;
l=1.8;
p.Training   = 'Multiscale';      % 'EIFEM'/'Multiscale', ('EIFisol')
p.Sampling   = 'Oversampling';  %'Isolated'/'Oversampling'
p.Inclusion  = 'Material';        % 'Material'/'Hole'
p.nelem      = 40;
meshName     = p.nelem+"x"+p.nelem;

%% DATA GENERATION

mR           = createReferenceMesh(p);
switch p.Training
    case 'Multiscale'
        p.Sampling      = 'Isolated';
        material  = createMaterial(mR,[1 1]);
        mesh            = mR;
        s.type          = 'discontinuous';
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
        [material,mT]    = createMaterial(mR,nS);
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
        data.material     = createMaterial(mR,[1 1]);
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

function mS=createStructuredMesh(p)
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

function uMesh=computeUnfittedMesh(bgMesh,levelSet)
    sUm.backgroundMesh = bgMesh;
    sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
    uMesh              = UnfittedMesh(sUm);
    uMesh.compute(levelSet);
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
            nS = [1 1 1]; %nx ny
            dI = [1 1 1];
        case 'Oversampling'
            nS = [3 3 3]; %nx ny
            dI = [2 2 2];
    end
end


function [material,m] = createMaterial(mesh,nSubdomains)
    mD = createMeshDomain(nSubdomains,mesh);
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

function mS = SmoothMesh(mR)
    s.coord    = mR.coord;
    s.connec   = mR.connec;

    xmin = min(s.coord(:,1));
    xmax = max(s.coord(:,1));
    ymin = min(s.coord(:,2));
    ymax = max(s.coord(:,2));
    zmin = min(s.coord(:,3));
    zmax = max(s.coord(:,3));
    tol=1e-10;

    leftNodes   = find(abs(s.coord(:,1)-xmin) < tol);
    rightNodes  = find(abs(s.coord(:,1)-xmax) < tol);
    bottomNodes = find(abs(s.coord(:,2)-ymin) < tol);
    topNodes    = find(abs(s.coord(:,2)-ymax) < tol);
    zNodes1     = find(abs(s.coord(:,3)-zmin) < tol);
    zNodes2     = find(abs(s.coord(:,3)-zmax) < tol);

    [~,iL] = sort(s.coord(leftNodes,2));
    leftNodes = leftNodes(iL);
    [~,iR] = sort(s.coord(rightNodes,2));
    rightNodes = rightNodes(iR);
    [~,iB] = sort(s.coord(bottomNodes,1));
    bottomNodes = bottomNodes(iB);
    [~,iT] = sort(s.coord(topNodes,1));
    topNodes = topNodes(iT);
    [~,iZ1] = sort(s.coord(zNodes1,3));
    zNodes1 = zNodes1(iZ1);
    [~,iZ2] = sort(s.coord(zNodes2,3));
    zNodes2 = zNodes2(iZ2);

    x = 0.5*(s.coord(bottomNodes,1)+s.coord(topNodes,1));
    s.coord(bottomNodes,1) = x;
    s.coord(topNodes,1)    = x;

    y = 0.5*(s.coord(leftNodes,2)+s.coord(rightNodes,2));
    s.coord(leftNodes,2)  = y;
    s.coord(rightNodes,2) = y;

    % z = 0.5*(s.coord(zNodes1,3)+s.coord(zNodes2,3));
    % s.coord(leftNodes,3)  = z;
    % s.coord(rightNodes,3) = z;

    mS = Mesh.create(s);
end