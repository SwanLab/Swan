%% DATA GENERATION

% This script has the purpose to create the necessary data related to the
% coarse training in order to obtain the NN and inputs for the
% preconditioner

clc; clear; close all;

%% INPUTS

 t1=0.03:0.05:0.499;
 t2=0.03:0.05:0.7;
%t1=0.1;
%t2=0.1;
p.Training   = 'EIFEM';      % 'EIFEM'/'Multiscale'
p.Sampling   = 'Oversampling';  %'Isolated'/'Oversampling'
p.nelem      = 50;
meshName     = p.nelem+"x"+p.nelem;


%% DATA GENERATION
for i=1:size(t1,2)
    for j = 1:size(t2,2)
        mR              = createReferenceMesh(p);
        g=computeLevelSet(t1(i),t2(j));
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
                s.dirichletFun  = obj.createDirichletFunction(bMesh);
                e  = ElasticHarmonicExtension(s);
                [T,lambda,K,Kcoarse] = e.solve();
                V   = material.computeVolume();
    
            case 'EIFEM'
                [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
                material     = createMaterial(mR,nS,'Material',g);
                s.mesh           = mR;
                s.material       = material;
                s.domainIndices  = dI;
                s.nSubdomains    = nS;            
                m= EIFEMTraining(s);
                data          = m.train();
                [data.material,mTr] = createMaterial(mR,[1 1],'Material',g);
                z = OfflineDataProcessor(data);
    
                EIFEoper = z.computeROMbasis();
                T        = EIFEoper.U;
                mesh     = data.mesh;
                Kcoarse  = EIFEoper.Kcoarse;
                vT      = mR.computeVolume();
                V        = Integrator.compute(mTr.designVariable.fun,mR,2);
                Vfr = V/vT;
        end
        tFrame  = t1(i);
        tCross  = t2(j);
    
        %Designa un nom per cada linea corresponent a un radi
        meshName=p.nelem+"x"+p.nelem;
        string = strrep("t1_"+num2str(t1(i), '%.2f'), ".", "_")+strrep("_t2_"+num2str(t2(j), '%.2f'), ".", "_")+"-"+meshName+".mat";
    
        % Guarda el .mat per cert radi
        FileName=fullfile('AbrilTFGfiles','Data2',p.Training,p.Sampling,string);
    
        switch p.Training
            case 'Multiscale'
                save(FileName,"T","Kcoarse","mesh","tFrame","tCross","V"); 
            case 'EIFEM'
                save(FileName, "EIFEoper","T","Kcoarse","mesh","tFrame","tCross","V","Vfr"); 
        end
    
    
    end
end


%% FUNCTIONS

function mS = createReferenceMesh(p)
    n =p.nelem;
    x1      = linspace(-1,1,n);
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
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymax,:)+[-delta,-0*delta];
    s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmax & s.coord(:,2)==obj.ymin,:)+[-delta,0*delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymax,:)+[delta,-0*delta];
    s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:) =...
        s.coord(s.coord(:,1)== obj.xmin & s.coord(:,2)==obj.ymin,:)+[delta,0*delta];

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
            nS = [1 1]; %nx ny
            dI = [1 1];
        case 'Oversampling'
            nS = [5 5]; %nx ny
            dI = [3 3];
    end
end

function g=computeLevelSet(t1,t2)
    gPar.type         = 'CrossedSquare';
    gPar.length       = 2;
    gPar.xCoorCenter  = 0;
    gPar.yCoorCenter  = 0;
    gPar.tFrame       = t1;
    gPar.tCross       = t2;
    g                 = GeometricalFunction(gPar);
end

function [material,m] = createMaterial(mesh,nSubdomains,inclusionType,g)
    s.mesh           = mesh;
    s.inclusionType  = inclusionType;
    s.nSubdomains    = nSubdomains;
    s.geomFun        = g;
    m = MaterialTraining(s);
    material = m.create();
end