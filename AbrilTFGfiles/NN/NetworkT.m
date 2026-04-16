% This script is the same one as the 1st version but it creates a loop and
% trains all the columns at once saving it all compactly in one variable to
% upload on the preconditioner.

clc;
clear;
close all;

%% Case parameters
p.Training  = 'Multiscale';   % 'EIFEM'/'Multiscale'
p.Sampling   ='Isolated';     %'Isolated'/'Oversampling'
p.Inclusion  ='Material';    %'Material'/'Hole'/'HoleRaul

%% Initialization of hyperparameters
pol_deg         = 5;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.32;
% hiddenLayers    = [224 250 280 300 280 250 224 200 150 100 72 50 20];
% hiddenLayers    = [300 350 400 450 450 400 350 300 200 100 72];
% hiddenLayers    = [100 150 200 250 200 150 100];
hiddenLayers    = [40 72 150 72];

%% INITIALIZATION 
% Store dataset file name
% s.fileName = fullfile('AbrilTFGfiles',"Data/Circle/",p.Training ,p.Inclusion,p.Sampling,'DataT.csv');
s.fileName = fullfile('AbrilTFGfiles','Data',"Cube",p.Training,'DataT.csv');

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 100000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';


% Select the model's features
s.xFeatures = 1:4;
s.yFeatures=[5:76];

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
    
% FileName=fullfile('AbrilTFGfiles',"Data/Circle/",p.Training,p.Inclusion,p.Sampling,"T_NN2.mat");
FileName=fullfile('AbrilTFGfiles',"Data",'Cube',p.Training,"T_NN.mat");
    save(FileName, "T_NN","pol_deg");
