%% T error Graphics Comparison

% This script has the purpose to visualize the error between the different
% options along the radius.

% 1. T NN 
% 2. HO FE aproximation
% 3. Hybrid SVD + NN

clc; clear;

%% LOAD DATA
p.Training  = 'EIFEM';            % 'EIFEM'/'Multiscale'
p.Inclusion  ='Material';         % 'Material'/'Hole'/'HoleRaul
p.Sampling   ='Isolated';         % 'Isolated'/'Oversampling'
p.nelem      = 20;
meshName    =  p.nelem+"x"+p.nelem;

% 1. NN
filePath = fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,"T_NN.mat");
load(filePath);
pol_deg1=pol_deg;

% 2. High Order function
HOname=fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,meshName,"HOfunction.mat");
load(HOname,"fT","deim");

% 3. SVD +NN
NNname=fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,meshName,"Q_NN.mat");
load(NNname);

U=deim.basis(:,1:10);
pol_deg2=pol_deg;


% Dataset
directory= fullfile("AbrilTFGfiles/Data",p.Training,p.Inclusion,p.Sampling,meshName);
files = dir(fullfile(directory, 'r0_*.mat'));
i=1;
for k = 1:1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    data=load(filePath);
    training.T(:,i) = data.T(:);  % This stores each file's contents in the cell array 'allData'
    training.mesh{i}=data.mesh;
    i=i+1;
end

%% TEST DATA GENERATION 

% test.r=0.025:0.05:0.999;
test.r=1e-6:0.025:0.951; 

test.T=zeros(size(training.T,1),size(test.r,2));
for i=1:size(test.r,2)
     mR=createReferenceMesh(p,test.r(i));
     test.mesh{i}=mR;
    switch p.Training
        case 'EIFEM'
            [nS,dI]      = defineNumberOfSubdomains(p.Sampling);
            material     = createMaterialTraining(mR, test.r(i),nS,p.Inclusion);
            s.mesh           = mR;
            s.r              = test.r(i);
            s.material       = material;
            s.domainIndices  = dI;
            s.nSubdomains    = nS;            
            m= EIFEMTraining(s);
            data          = m.train();
            data.material = createMaterialTraining(mR, test.r(i) ,[1 1],p.Inclusion);
            z = OfflineDataProcessor(data);
            EIFEoper = z.computeROMbasis();
            test.T(:,i)=EIFEoper.U(:);
        case 'Multiscale'
            p.Sampling = 'Isolated';
            material   = createMaterialTraining(mR, test.r(i),[1 1],p.Inclusion);
            mesh       = mR;
            bMesh      = mesh.createSingleBoundaryMesh();
            s.mesh=bMesh;
            cf=CoarseFunctions(s);
            f = cf.getAnalytical();
            s.mesh          = mesh;
            s.uFun          = LagrangianFunction.create(mesh, mesh.ndim, 'P1');
            s.lambdaFun     = LagrangianFunction.create(bMesh,mesh.ndim, 'P1');
            s.material      = material;
            s.dirichletFun  = f;
            e  = ElasticHarmonicExtension(s);
            [T,~,~,~] = e.solve();
            test.T(:,i)= T(:);
    end
end


%% RECONSTRUCT T DATA
training.r=0:0.05:0.999;
mesh=training.mesh{1}; % Totes tenen mateix num de dofs;

training.T1= zeros(mesh.nnodes*mesh.ndim*8,length(training.r));
training.T2= zeros(mesh.nnodes*mesh.ndim*8,length(training.r));
training.T3= zeros(mesh.nnodes*mesh.ndim*8,length(training.r));

test.T1= zeros(mesh.nnodes*mesh.ndim*8,length(test.r));
test.T2= zeros(mesh.nnodes*mesh.ndim*8,length(test.r));
test.T3= zeros(mesh.nnodes*mesh.ndim*8,length(test.r));

% 1. NN
for i=1:length(training.r)
    aux=computeT_NN(training.mesh{i},training.r(i),T_NN,pol_deg1);
    training.T1(:,i)=aux(:);
end

for i=1:length(test.r)
    aux=computeT_NN(test.mesh{i},test.r(i),T_NN,pol_deg1);
    test.T1(:,i)=aux(:);
end


% 2. High Order function
for i=1:length(training.r)
    aux=fT(training.r(i));
    training.T2(:,i)=aux(:);
end

for i=1:length(test.r)
    aux=fT(test.r(i));
    test.T2(:,i)=aux(:);
end


% 3. SVD + NN
for i=1:length(training.r)
    aux=computeT_Hybrid(basis,training.r(i),Q_NN,pol_deg2);
    training.T3(:,i)=aux(:);
end

for i=1:length(test.r)
    aux=computeT_Hybrid(basis,test.r(i),Q_NN,pol_deg2);
    test.T3(:,i)=aux(:);
end


%% CALCULATE ERROR and test

training.err1=vecnorm(abs(training.T-training.T1))/norm(training.T);
training.err2=vecnorm(abs(training.T-training.T2))/norm(training.T);
training.err3=vecnorm(abs(training.T-training.T3))/norm(training.T);

test.err1=vecnorm(abs(test.T-test.T1))/norm(test.T);
test.err2=vecnorm(abs(test.T-test.T2))/norm(test.T);
test.err3=vecnorm(abs(test.T-test.T3))/norm(test.T);

%% PLOT ERROR

figure
plot(training.r,training.err1,training.r,training.err2,training.r,training.err3,LineWidth=1.5);
legend("NN","HO","SVD+NN");
title("Training Error vs r");
xlabel('r');
ylabel('error');


%% PLOT TEST

figure
plot(test.r,test.err1,test.r,test.err2,test.r,test.err3,LineWidth=1.5);
legend("NN","HO","SVD+NN");
title("Test Error vs r ");
xlabel('r');
ylabel('error');

figure
plot(test.r,test.err2,LineWidth=1.5);
legend("SVD+HO")
title("Test Error vs r ");
xlabel('r');
ylabel('error');


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
