clc;
clear;
close all;

%% Case parameters
p.Training  = 'EIFEM';        % 'EIFEM'/'Multiscale'
p.Inclusion  ='Material';    %'Material'/'Hole'/'HoleRaul
p.Sampling   ='Isolated';     %'Isolated'/'Oversampling'

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.1;
hiddenLayers    = [10 20 40 80 100 80 40];
% hiddenLayers    = [24 30 40 80 100 80 40];

%% INITIALIZATION 
% Store dataset file name
% s.fileName = fullfile('AbrilTFGfiles',"Data",p.Training,p.Inclusion,p.Sampling,'DataK.csv');
s.fileName = fullfile('AbrilTFGfiles',"Data",p.Training ,'Lattice','DataK.csv');

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 1500000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = [1:2];
s.yFeatures = [3:1:38];
cHomogIdxs = [11, 12, 22, 33];

% Load data
%data   = cHomogData(s);
data   = Data(s);
s.data = data;

% Train the model
K_NN = OptimizationProblemNN(s);
K_NN.solve();
K_NN.plotCostFnc();
MSETrain    = immse(K_NN.computeOutputValues(data.Xtrain), data.Ytrain);

string ="K_NN.mat";
% FileName=fullfile('AbrilTFGfiles',"Data",p.Training,p.Inclusion,p.Sampling,string);
%     save(FileName, "K_NN");

FileName=fullfile('AbrilTFGfiles',"Data",p.Training,'Lattice',"K_NN.mat");
save(FileName, "K_NN","pol_deg");

%%% Plot surface
%
%% Load dataset from specified path
%filePath = fullfile('AbrilTFGfiles', s.fileName);
%tempData = readmatrix(filePath);
%
%real = tempData(:,s.yFeatures);
%predicted = zeros(size(real));
%
%for i = 1:size(real,1)
%    predicted(i,:) = K_NN.computeOutputValues(tempData(i,s.xFeatures));
%end
%
%difference = real-predicted;
%disp(norm(difference));
