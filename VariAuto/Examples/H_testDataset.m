%% Testing with a sample nn
clc;
clear;
close all;
%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 0;
lambda          = 0;
learningRate    = 0.1;
momentum        = 0.9;
batch           = 200;
hiddenlayers    = [];

%% Loading of files/datasets
fileN = 'regression_dataset.csv';
s.features = 1:1;
s.fileName        = fileN;
s.testRatio       = testratio;
s.polynomialOrder = pol_deg;
data  = Data(s);

% data.Ytest = data.Ytest(~isnan(data.Ytest(:,1)),:);
% data.Xtest = data.Xtest(~isnan(data.Ytest(:,1)),:);
% data.Xtrain = data.Xtrain(~isnan(data.Ytrain(:,1)),:);
% data.Ytrain = data.Ytrain(~isnan(data.Ytrain(:,1)),:);
% data.Ntest = size(data.Ytest,1);

%% Create Network and trainer Objects
structure = [data.nFeatures,hiddenlayers,data.nLabels];

%% Run Optimization Problem
p.data            = data;
p.structure       = structure;
p.optimizerParams.learningRate = learningRate;
p.costParams.lambda = lambda;
p.networkParams.hiddenLayers = hiddenlayers;
p.networkParams.costType     = 'L2';
p.networkParams.HUtype       = 'None';
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











