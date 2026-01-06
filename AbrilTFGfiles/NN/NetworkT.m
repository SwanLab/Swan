% This script is the same one as the 1st version but it creates a loop and
% trains all the columns at once saving it all compactly in one variable to
% upload on the preconditioner.

clc;
clear;
close all;

%% Case parameters
p.Training  = 'EIFEM';        % 'EIFEM'/'Multiscale'
p.Sampling   ='Oversampling';     %'Isolated'/'Oversampling'
p.Inclusion  ='Material';    %'Material'/'Hole'/'HoleRaul

%% Initialization of hyperparameters
pol_deg         = 9;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.01;
hiddenLayers    = [224 250 280 300 280 250 224 200 150 100 72 50 20];


%% INITIALIZATION 
% Store dataset file name
s.fileName = fullfile('AbrilTFGfiles',"Data",p.Training ,p.Inclusion,p.Sampling,'DataT.csv');

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 10000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';


% Select the model's features
s.xFeatures = 1:3;
s.yFeatures=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19];


%% Initialization of variables to save
MSETrain=zeros(1,8);
comparison=cell(1,8);

%% Loop for the 8 coarse modes

% Load data
data   = Data(s);
s.data = data;

% Train the model
T_NN = OptimizationProblemNN(s);
T_NN.solve();
T_NN.plotCostFnc();
    
FileName=fullfile('AbrilTFGfiles',"Data",p.Training,p.Inclusion,p.Sampling,"T_NN.mat")
    save(FileName, "T_NN","pol_deg");
