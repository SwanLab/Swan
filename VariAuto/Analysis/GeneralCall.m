%% Testing with a sample nn
clc;
clear;
close all;
addpath ../Codes;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;  
lambda          = 0;
learningRate    = 0.1;
momentum        = 0.9;
batch           = 200;
hiddenlayers    = [10,15];

%% Loading of files/datasets
datasets = load("../Codes/datasets.mat").datasets1;
disp('Datsets available:')
for i = 1:length(datasets)
    fprintf('%d - %s \n',i,datasets(i))
end
fileN = datasets(input('Choose: '));
data  = Data(fileN,testratio,pol_deg);

%% Create Network and trainer Objects
structure = [data.nFeatures,hiddenlayers,data.nLabels];
% network   = Network(data,structure);
% network = Network(data,structure,'-loglikelihood','ReLU','softmax',lambda);

%% Run Optimization Problem
optProblem   = optimizationProblem(data,structure,learningRate);
% opt.optTolerance  = 1*10^-8; opt.maxevals      = 100;
% opt.maxepochs     = 100    ; opt.earlyStop     = 10;
% opt.time          = Inf([1,1]); opt.fv         = 10^-4;
% nplt              = 1;
% optimizer       = Trainer.create(network,'SGD',learningRate,momentum,batch,opt,'static',nplt);

%% RUN & Possible functions
data.plotCorrMatrix();
% network.plotBoundary('contour'); Amb errors dins?
optProblem.plotConections();
optProblem.plotConfusionMatrix();


