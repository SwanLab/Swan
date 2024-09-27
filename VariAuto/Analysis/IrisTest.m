%% Testing with IRIS DATASET
clc;
clear;
close all;
%addpath ../Codes;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;  
lambda          = 0;
learningRate    = 0.1;
hiddenLayers    = [3,2,1,2,3];

%% INITIALIZATION 
% try different feature combination, pairs of features enable the possibility of plotting boundaries
s.fileName = 'Iris.csv';
s.polynomialOrder = 1;
s.testRatio       = 30;
data = Data(s);

s.networkParams.hiddenLayers    = hiddenLayers;
s.data                          = data;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;


opt = OptimizationProblem(s);
opt.solve();
opt.plotConfusionMatrix();
%opt.plotSurface();

if data.nFeatures == 2  %If you want to be asked for Features change it in "Data" Class
    opt.plotBoundary('contour');
end
