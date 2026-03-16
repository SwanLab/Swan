%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

% r=0:0.05:0.999;
r=0.5;

p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale'
p.Sampling   = 'Oversampling';  %'Isolated'/'Oversampling'
p.nelem      = 20;
meshName     = p.nelem+"x"+p.nelem;

%% DATA GENERATION
for j = 1:size(r,2)
    radius = r(j);
        mR              = createReferenceMesh(p);
        g=computeLevelSet(radius);
        switch p.Training
            case 'Multiscale'
                p.Sampling = 'Isolated';
                material   = createMaterial(mR,[1 1],'Material',g);
                mesh       = mR;
                bMesh      = mesh.createSingleBoundaryMesh();
                s.mesh          = mesh;
                s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
                s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
                s.material      = material;
                s.dirichletFun  = createDirichletFunction(bMesh);
                e  = ElasticHarmonicExtension(s);
                [T,lambda,K,Kcoarse] = e.solve();
                V  = Integrator.compute(mTr.designVariable.fun,mR,2);
    
            case 'EIFEM'
                [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
                [material,mT]     = createMaterial(mR,nS,'Material',g);
                s.mesh           = mR;
                s.material       = material;
                s.domainIndices  = dI;
                s.nSubdomains    = nS;  
                s.levelSet = mT.levelSet;
                s.unfittedMesh = mT.unfittedMesh;
                m= EIFEMTraining(s);
                data          = m.train();
                [data.material,mTr] = createMaterial(mR,[1 1],'Material',g);
                z = OfflineDataProcessor(data);
    
                EIFEoper = z.computeROMbasis();
                T        = EIFEoper.U;
                mesh     = data.mesh;
                Kcoarse  = EIFEoper.Kcoarse;
                vT       = mR.computeVolume();
                V        = Integrator.compute(mTr.designVariable.fun,mR,2);
                Vfr = V/vT;
        end
    
        R        = r(j);
        %Designa un nom per cada linea corresponent a un radi
        meshName=p.nelem+"x"+p.nelem;
        string = strrep("r"+num2str(r(j), '%.4f'), ".", "_")+"-"+meshName+".mat";
    
        % Guarda el .mat per cert radi
        FileName=fullfile('AbrilTFGfiles','Data2',p.Training,"Sphere",string);
    
        switch p.Training
            case 'Multiscale'
                save(FileName,"T","Kcoarse","mesh","R"); 
            case 'EIFEM'
                save(FileName, "EIFEoper","T","Kcoarse","mesh","R"); 
        end
    
end


%% FUNCTIONS

function mS = createReferenceMesh(p)
    n =p.nelem;
    m = TetraMesh(1,1,1,n,n,n);

    s.coord=m.coord;
    s.connec=m.connec;
    
    maxC= max(s.coord);
    minC = min(s.coord);
    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:) =s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];
    
    s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== maxC(1) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];
    
    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];
    
    s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,1)== minC(1) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];
    
    s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];
    
    s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,2)== maxC(2) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];
    
    s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:) =...
        s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==maxC(3),:)-[0,0,1e-5];
    
    s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:) =...
        s.coord(s.coord(:,2)== minC(2) & s.coord(:,3)==minC(3),:)+[0,0,1e-5];
    
    mS = Mesh.create(s);

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

function g=computeLevelSet(r)
    gPar.type         = 'Sphere';
    gPar.radius       = r;
    gPar.xCoorCenter  = 0;
    gPar.yCoorCenter  = 0;
    gPar.zCoorCenter  = 0;
    g                 = GeometricalFunction(gPar);
end

function [material,m] = createMaterial(mesh,nSubdomains,inclusionType,g)
    s.mesh       = mesh;
    s.inclusionType  = inclusionType;
    s.nSubdomains    = nSubdomains;
    s.geomFun       = g;
    m = MaterialTraining(s);
    material = m.create();
end