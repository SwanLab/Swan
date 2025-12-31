%SVD FOR THE NN OPTION 3

clc;
clear;
close all;

%% Load Data
p.nelem     =  10;
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

fileName=fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,meshName,"SVD.mat");
save(fileName,"U","S","V","r","basis","nBasis");
QFileName = fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,meshName,"DataQ.csv");
writematrix(table,QFileName);
