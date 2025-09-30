clc;
clear;
close all;
%% Initialization of hyperparameters
%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 50;
lambda          = 0;
learningRate    = 0.01;

hiddenLayers    = 256 .* ones(1,16);
%% INITIALIZATION 
% Store dataset file name
s.fileName = 'AlbertTFG files/comparisons/k case 3/K_data case 3 100x100.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 10^9; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'tanh';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = 1;
s.yFeatures = 2:1:65;

% Load data
%data   = cHomogData(s);
data   = Data(s);
s.data = data;

% Train the model
opt = OptimizationProblemNN(s);
opt.solve();
opt.plotCostFnc();
MSETrain    = immse(opt.computeOutputValues(data.Xtrain), data.Ytrain);

% %% Plot surface
% 
% % Load dataset from specified path
% filePath = fullfile('AlbertTFG files/mat files/Full/csvFiles_Case3', s.fileName);
% tempData = readmatrix(filePath);
% 
% real = tempData(:, s.yFeatures);
% predicted = zeros(size(real));
% 
% for i = 1:size(real,1)
%     predicted(i, :) = opt.computeOutputValues(tempData(i, s.xFeatures));
% end
% 
% difference = real-predicted;

% figure
% bar3(abs(difference))
% title("Abs of real vs predicted")
% xlabel("Position i in array")
% ylabel("Position j in array")
% 
% 
% predictedHalf = zeros(size(real));
% 
% for i = 1:size(real,1)
%     predictedHalf(i, :) = opt.computeOutputValues(0.5.*tempData(i, s.xFeatures));
% end
% differenceHalf = 0.5.*real-predictedHalf;
% 
% figure
% bar3(abs(differenceHalf))
% title("Abs of real vs predicted when input is 0.5")
% xlabel("Position i in array")
% ylabel("Position j in array")
