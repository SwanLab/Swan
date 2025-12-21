%SVD FOR THE NN OPTION 3

clc;
clear;
close all;

%% Load Data
p.nelem     =  20;
p.Inclusion = 'HoleRaul';         % 'Hole'/'Material'/'HoleRaul'
p.Sampling  = 'Isolated';     % 'Isolated'/'Oversampling'
meshName    = p.nelem+"x"+p.nelem;

files = dir(fullfile("AbrilTFGfiles/Data/",p.Inclusion,p.Sampling,meshName, 'r0_*.mat'));

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

fileName=fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,"SVD.mat");
save(fileName,"U","S","V","r");
QFileName = fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,"DataQ.csv");
writematrix(table,QFileName);


%% Export NN paraview file
NNname=fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,"Q_NN.mat");
load(NNname);


rad=0.3;
rFull = Data.buildModel(rad,6);
qNN=Q_NN.computeOutputValues(rFull).';


T_NN=basis*qNN;
u=reshape(T_NN,[],8);

exportT_weakInclusion(u,rad, data.mesh,"TestNN");