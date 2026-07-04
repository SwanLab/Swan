%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale'
p.Sampling   = 'Isolated';  %'Isolated'/'Oversampling'
p.Inclusion  = 'Hole';        % 'Material'/'Hole'/'MeshRaul'
p.nelem      = 40;

%% DATA GENERATION
mR              = createReferenceMesh(p);
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
        s.type           = 'continuous';
        s.Coarseorder    = 2;
        m= EIFEMTraining(s);
        data             = m.train();
        [data.material]  = createMaterial(mR,[1 1]);
        data.dirac       = true;
        data.type        = s.type;
        data.Coarseorder = 2;
        z = OfflineDataProcessor(data);

        EIFEoper = z.computeROMbasis();
        T        = EIFEoper.U;
        mesh     = data.mesh;
        Kcoarse  = EIFEoper.Kcoarse;
end

string = "AuxeticIsolatedClassic.mat";

% Guarda el .mat per cert radi
FileName=fullfile('AbrilTFGfiles','Data',"Auxetic",string);

switch p.Training
    case 'Multiscale'
        save(FileName,"T","Kcoarse");
    case 'EIFEM'
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
    
        case 'MeshRaul'
            % filename = 'DEF_Q4auxL_1.mat';
            filename = 'DEF_Q4porL_1.mat';
            load(filename,'EIFEoper');
            s.coord    = EIFEoper.MESH.COOR;
            s.connec   = EIFEoper.MESH.CN;
    
            obj.xmin = min(s.coord(:,1));
            obj.xmax = max(s.coord(:,1));
            obj.ymin = min(s.coord(:,2));
            obj.ymax = max(s.coord(:,2));
    
            delta = 0;
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
end

function mS = SmoothMesh(mR)
    s.coord    = mR.coord;
    s.connec   = mR.connec;

    xmin = min(s.coord(:,1));
    xmax = max(s.coord(:,1));
    ymin = min(s.coord(:,2));
    ymax = max(s.coord(:,2));

    tol=1e-10;

    leftNodes   = find(abs(s.coord(:,1)-xmin) < tol);
    rightNodes  = find(abs(s.coord(:,1)-xmax) < tol);
    bottomNodes = find(abs(s.coord(:,2)-ymin) < tol);
    topNodes    = find(abs(s.coord(:,2)-ymax) < tol);

    [~,iL] = sort(s.coord(leftNodes,2));
    leftNodes = leftNodes(iL);
    [~,iR] = sort(s.coord(rightNodes,2));
    rightNodes = rightNodes(iR);
    [~,iB] = sort(s.coord(bottomNodes,1));
    bottomNodes = bottomNodes(iB);
    [~,iT] = sort(s.coord(topNodes,1));
    topNodes = topNodes(iT);
    y = 0.5*(s.coord(leftNodes,2)+s.coord(rightNodes,2));
    s.coord(leftNodes,2)  = y;
    s.coord(rightNodes,2) = y;

    x = 0.5*(s.coord(bottomNodes,1)+s.coord(topNodes,1));
    s.coord(bottomNodes,1) = x;
    s.coord(topNodes,1)    = x;
    mS = Mesh.create(s);


end

function dF = createDirichletFunction(bMesh,type)
s.mesh = bMesh;
s.order= 2;
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
            nS = [3 3]; %nx ny
            dI = [2 2];
    end
end


function [material] = createMaterial(mesh,nS)
    mD      = createMeshDomain(mesh,nS);
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

function mS = createStructuredMesh(p)
    n =p.nelem;
    x1      = linspace(-1.5,1.5);
    x2      = linspace(-1,1,n);
    [xv,yv] = meshgrid(x1,x2);
    [F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
    s.coord  = V(:,1:2);
    s.connec = F;
    obj.xmin = min(x1);            
    obj.xmax = max(x1);
    obj.ymin = min(x2);
    obj.ymax = max(x2);
    delta = 1e-9;
    s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:) =...
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-delta];
    s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,+delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[+delta,-delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[+delta,+delta];
    mS = Mesh.create(s);
end

function levelSet = createLevelSetFunction(bgMesh)
    gPar.type           = 'Auxetic';
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