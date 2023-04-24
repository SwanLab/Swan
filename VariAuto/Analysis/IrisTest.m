%% Testing with IRIS DATASET
clc;
clear;
close all;
addpath ../Codes;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;  
lambda          = 0;
learningRate    = 0.1;
hiddenlayers    = [2,3];

%% INITIALIZATION 
% try different feature combination, pairs of features enable the possibility of plotting boundaries
data      = Data('../Datasets/Iris.csv',30,1);
structure = [data.nFeatures,hiddenlayers,data.nLabels];
optProblem   = optimizationProblem(data,structure,learningRate,lambda);

optProblem.plotConfusionMatrix();

if data.nFeatures == 2  %If you want to be asked for Features change it in "Data" Class
    optProblem.plotBoundary('contour');
end
