%% T error Graphics Comparison

% This script has the purpose to visualize the error between the different
% options along the radius.

% 1. T NN 
% 2. HO FE aproximation
% 3. Hybrid SVD + NN

clc; clear;

%% LOAD DATA
p.Training  = 'EIFEM';            % 'EIFEM'/'Multiscale'
p.Sampling   ='Oversampling';     %'Isolated'/'Oversampling'
p.Inclusion  ='HoleRaul';         %'Material'/'Hole'/'HoleRaul
p.nelem      = 20;
meshName    =  p.nelem+"x"+p.nelem;

% 1. NN
filePath = fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,"T_NN.mat");
% filePath = fullfile("AbrilTFGfiles","Data","Multiscale","T_NN.mat");
load(filePath);
pol_deg1=pol_deg;

% 2. High Order function
HOname=fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,meshName,"HOfunction.mat");
% HOname=fullfile("AbrilTFGfiles","Data","Multiscale","HOfunction.mat");
load(HOname,"fT","deim");

% 3. SVD +NN
NNname=fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,meshName,"Q_NN.mat");
% NNname=fullfile("AbrilTFGfiles","Data","Multiscale","Q_NN.mat");
load(NNname);

U=deim.basis(:,1:10);
pol_deg2=pol_deg;


% Dataset
directory= fullfile("AbrilTFGfiles/Data",p.Training,p.Inclusion,p.Sampling,meshName);
% directory= fullfile("AbrilTFGfiles/Data/Multiscale/");
files = dir(fullfile(directory, 'r0_*.mat'));
i=1;
for k = 1:1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    data=load(filePath);
    training.T(:,i) = data.T(:);  % This stores each file's contents in the cell array 'allData'
    i=i+1;
end

mesh=data.mesh;

%% TEST DATA GENERATION 

test.r=0.025:0.05:0.999;

test.T=zeros(size(training.T,1),size(test.r,2));
for i=1:size(test.r,2)
    % EIFEM:
    %data=OversamplingTraining(mesh,test.r(i),p);
    %z=OfflineDataProcessor(data);
    %EIFEoper = z.computeROMbasis();
    %test.T(:,i)=EIFEoper.U(:);
    
    % MULTISCALE:
    [~, ~, T, ~, ~,~] = IsolatedTraining(test.r(i),p.nelem);
    test.T(:,i)= T(:);
end


%% RECONSTRUCT T DATA
training.r=0:0.05:0.999;

training.T1= zeros(mesh.nnodes*mesh.ndim*8,length(training.r));
training.T2= zeros(mesh.nnodes*mesh.ndim*8,length(training.r));
training.T3= zeros(mesh.nnodes*mesh.ndim*8,length(training.r));

test.T1= zeros(mesh.nnodes*mesh.ndim*8,length(test.r));
test.T2= zeros(mesh.nnodes*mesh.ndim*8,length(test.r));
test.T3= zeros(mesh.nnodes*mesh.ndim*8,length(test.r));

% 1. NN
for i=1:length(training.r)
    aux=computeT_NN(mesh,training.r(i),T_NN,pol_deg1);
    training.T1(:,i)=aux(:);
end

for i=1:length(test.r)
    aux=computeT_NN(mesh,test.r(i),T_NN,pol_deg1);
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
