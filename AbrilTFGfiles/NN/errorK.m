%% T error Graphics Comparison

% This script has the purpose to visualize the error between the different
% options along the radius.

% 1. K NN 
% 2. HO FE aproximation

clc; clear;

%% LOAD DATA
p.Training  = 'EIFEM';            % 'EIFEM'/'Multiscale'
p.Sampling   ='Oversampling';         %'Isolated'/'Oversampling'
p.Inclusion  ='Material';         %'Material'/'Hole'/'HoleRaul
p.nelem      = 20;
meshName    =  p.nelem+"x"+p.nelem;

% 1. NN
filePath = fullfile("AbrilTFGfiles","Data",'Circle/',p.Training,p.Inclusion,p.Sampling,"K_NN.mat");
load(filePath);

% 2. High Order function
HOname=fullfile("AbrilTFGfiles","Data","Circle/",p.Training,p.Inclusion,p.Sampling,meshName,"parametrizedEIFEM.mat");
load(HOname,"EIFEoper");
fK = EIFEoper.Kcoarse;

% Dataset
directory= fullfile("AbrilTFGfiles/Data","Circle/",p.Training,p.Inclusion,p.Sampling,meshName);
files = dir(fullfile(directory, 'r0_*.mat'));
i=1;
for k = 1:1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    data=load(filePath);
    training.Kcoarse(:,i) = data.Kcoarse(:);  % This stores each file's contents in the cell array 'allData'
    training.mesh{i}=data.mesh;
    i=i+1;
end

%% TEST DATA GENERATION 

% test.r=0.025:0.05:0.999;
test.r=1e-6:0.025:0.951; 
test.Kcoarse=zeros(size(training.Kcoarse,1),size(test.r,2));

for i=1:size(test.r,2)
     mR=createReferenceMesh(p,test.r(i));
     test.mesh{i}=mR;
     g=computeLevelSet(test.r(i));
    switch p.Training
        case 'EIFEM'
            [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
            material     = createMaterial(mR,nS,p.Inclusion,g);
            s.mesh           = mR;
            s.r              = test.r(i);
            s.material       = material;
            s.domainIndices  = dI;
            s.nSubdomains    = nS;            
            m= EIFEMTraining(s);
            data          = m.train();
            data.material = createMaterial(mR,[1 1],p.Inclusion,g);
            z = OfflineDataProcessor(data);
            EIFEoper = z.computeROMbasis();
            test.Kcoarse(:,i)=EIFEoper.Kcoarse(:);
        case 'Multiscale'
            p.Sampling = 'Isolated';
            material   = createMaterial(mR,[1 1],p.Inclusion,g);
            mesh       = mR;
            bMesh      = mesh.createSingleBoundaryMesh();
            s.mesh=bMesh;
            s.type='continuous';
            cf=CoarseFunctions(s);
            f = cf.getAnalytical();
            s.mesh          = mesh;
            s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
            s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
            s.material      = material;
            s.dirichletFun  = f;
            e  = ElasticHarmonicExtension(s);
            [~,~,~,Kc] = e.solve();
            test.Kcoarse(:,i)= Kc(:);
    end
end


%% RECONSTRUCT T DATA
training.r=0:0.05:0.999;
mesh=training.mesh{1}; % Totes tenen mateix num de dofs;

training.K1= zeros(8*8,length(training.r));
training.K2= zeros(8*8,length(training.r));

test.T1= zeros(8*8,length(test.r));
test.T2= zeros(8*8,length(test.r));

% 1. NN
for i=1:length(training.r)
    aux=computeKcoarse_NN(K_NN,training.r(i));
    training.K1(:,i)=aux(:);
end

for i=1:length(test.r)
    aux=computeKcoarse_NN(K_NN,test.r(i));
    test.T1(:,i)=aux(:);
end


% 2. High Order function
for i=1:length(training.r)
    aux=fK(training.r(i));
    training.K2(:,i)=aux(:);
end

for i=1:length(test.r)
    aux=fK(test.r(i));
    test.T2(:,i)=aux(:);
end


%% CALCULATE ERROR and test

training.err1=vecnorm(abs(training.Kcoarse-training.K1))/norm(training.Kcoarse);
training.err2=vecnorm(abs(training.Kcoarse-training.K2))/norm(training.Kcoarse);

test.err1=vecnorm(abs(test.Kcoarse-test.T1))/norm(test.Kcoarse);
test.err2=vecnorm(abs(test.Kcoarse-test.T2))/norm(test.Kcoarse);

%% PLOT ERROR

figure
tiledlayout(1,2)
pos= [572,507,978,371];
set(gcf, 'Position', pos) 
nexttile
plot(training.r,training.err1,LineWidth=1.5);
% legend("NN","HO");
title("NN training Error vs r");
xlabel('r');
ylabel('error');
ylim ([0 0.03]);

%% PLOT TEST

nexttile
plot(test.r,test.err1,LineWidth=1.5);
% legend("NN","HO");
title("NN test Error vs r ");
xlabel('r');
ylabel('error');
ylim ([0 0.03]);

% %% Prova plots
% 
% 
% ri = linspace(min(test.r), max(test.r), 500);   % más puntos
% err2i = interp1(test.r, test.err2, ri, 'spline');     % o 'pchip'
% 
% figure
% plot(test.r, test.err2, 'o', ri, err2i, '-')
% 
% err2_s = smoothdata(y, 'sgolay');


%% FUNCTIONS

function mS = createReferenceMesh(p,r)
    
    switch p.Inclusion
        case 'Material'
            mS = createStructuredMesh(p);
        case 'Hole'
            mS = createStructuredMesh(p);
            lvSet    = createLevelSetFunction(mS,r);
            uMesh    = computeUnfittedMesh(mS,lvSet);
            mS       = uMesh.createInnerMesh();
        case 'HoleRaul'
            switch p.nelem
                case 10
                    mS=mesh_rectangle_via_triangles(r,1,-1,1,-1,7,6,0,0);   % 10x10
                case 20
                    mS=mesh_rectangle_via_triangles(r,1,-1,1,-1,15,12,0,0); % 20x20
                case 50
                    mS=mesh_rectangle_via_triangles(r,1,-1,1,-1,34,35,0,0)  % 50x50
            end
    end
end


function mS = createStructuredMesh(p)
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

function levelSet = createLevelSetFunction(bgMesh,r)
    sLS.type        = 'CircleInclusion';
    sLS.xCoorCenter = 0;
    sLS.yCoorCenter = 0;
    sLS.radius      = r;
    g               = GeometricalFunction(sLS);
    lsFun           = g.computeLevelSetFunction(bgMesh);
    levelSet        = lsFun.fValues;
end

function uMesh = computeUnfittedMesh(bgMesh,levelSet)
    sUm.backgroundMesh = bgMesh;
    sUm.boundaryMesh   = bgMesh.createBoundaryMesh();
    uMesh              = UnfittedMesh(sUm);
    uMesh.compute(levelSet);
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

function material = createMaterial(mesh,nSubdomains,inclusionType,g)
    s.mesh           = mesh;
    s.inclusionType  = inclusionType;
    s.nSubdomains    = nSubdomains;
    s.geomFun        = g;
    m = MaterialTraining(s);
    material = m.create();
end

function g=computeLevelSet(r)
    gPar.type         = 'Circle';
    gPar.radius       = r;
    gPar.xCoorCenter  = 0;
    gPar.yCoorCenter  = 0;
    g                 = GeometricalFunction(gPar);
end