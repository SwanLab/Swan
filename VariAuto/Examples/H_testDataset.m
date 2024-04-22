%% Testing with a sample nn
clc;
clear;
close all;
%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;  
lambda          = 0.1;
learningRate    = 0.1;
momentum        = 0.9;
batch           = 200;
hiddenlayers    = [10,15];

%% Loading of files/datasets
fileN = 'testDataset.csv';
s.features = 1:3;
s.fileName        = fileN;
s.testRatio       = testratio;
s.polynomialOrder = pol_deg;
data  = Data(s);

%% Create Network and trainer Objects
structure = [data.nFeatures,hiddenlayers,data.nLabels];

%% Run Optimization Problem
p.data            = data;
p.structure       = structure;
p.optimizerParams.learningRate = learningRate;
p.costParams.lambda = lambda;
p.networkParams.hiddenLayers = hiddenlayers;
p.networkParams.costType     = 'L2';
p.networkParams.HUtype       = 'ReLU';
p.networkParams.OUtype       = 'None';
optProblem   = OptimizationProblem(p);

optProblem.solve();

%% Additional plots
optProblem.plotCostFnc();

%% Error
% Evaluation for the train data
errTrain = optProblem.computeError(data.Xtrain, data.Ytrain);
% Evaluation for the test data
errTest = optProblem.computeError(data.Xtest, data.Ytest);











