%SVD FOR THE NN OPTION 3

clc;
clear;
close all;

%% Load Data
directory = './AbrilTFGfiles/Data/50x50';
files     = dir(fullfile(directory, 'UL_*.mat'));

r=zeros(length(files),1);
for k = 1:1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    data=load(filePath);
    if k==1
        T=zeros(data.mesh.nnodes*data.mesh.ndim*8,length(files));
    end
    T(:,k) = data.T(:);
    r(k,1)=data.R;
end

%% Data to train
[U,S,V]    = svd(T,"econ");
nBasis=10;

basis = U(:,1:nBasis);
sValues = S(1:nBasis,1:nBasis);
q = V(:,1:nBasis)*sValues;

table=[r q];

%fileName=fullfile("AbrilTFGfiles","SVD","SVD.mat");
%save(fileName,"U","S","V","r");
%QFileName = fullfile('AbrilTFGfiles', 'DataQ.csv');
%writematrix(table,QFileName);


%% Export NN paraview file
NNname=fullfile("AbrilTFGfiles/NN/Q_NN2.mat");
load(NNname);


rad=0.3;
rFull = Data.buildModel(rad,6);
qNN=Q_NN.computeOutputValues(rFull).';


T_NN=basis*qNN;

u=reshape(T_NN,[],8);

% Export vtu file
 z.mesh      = data.mesh;
 z.order     = 'P1';
 
 for i=1:8
   z.fValues   = reshape(u(:,i),[data.mesh.ndim,data.mesh.nnodes])';
   uFeFun = LagrangianFunction(z);%
   fileName = strrep("r" + num2str(rad), '.', '_')+ "_SVD_NN_Training" +num2str(i);
   centroids=computeCentroid(data.mesh);
   CoarsePlotSolution(uFeFun, data.mesh,[],fileName, rad, centroids);
 end

function centroids=computeCentroid(mesh)
    x0=mean(mesh.coord(:,1));
    y0=mean(mesh.coord(:,2));
    centroids = [x0,y0];
end
    