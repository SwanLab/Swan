clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 5;
lambda          = 0.0;
learningRate    = 0.2;
hiddenLayers    = 16 .* ones(1, 1);

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'AlbertTFG files/comparisons/k case 1/NewL/K_r0.9000-20x20.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 5000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = 1:1:8;
s.yFeatures = 9:1:16;

% Load data
%data   = cHomogData(s);
data   = Data(s);
s.data = data;

% Train the model
opt = OptimizationProblemNN(s);
opt.solve();
opt.plotCostFnc();
MSETrain    = immse(opt.computeOutputValues(data.Xtrain), data.Ytrain);

%% Plot surface

% Load dataset from specified path
filePath = fullfile(s.fileName);
tempData = readmatrix(filePath);

real = tempData(:, s.yFeatures);
predicted = zeros(size(real));

for i = 1:size(real,1)
    predicted(i, :) = opt.computeOutputValues(tempData(i, s.xFeatures));
end

difference = real-predicted;

figure
bar3c(abs(difference));
zlabel("abs(diff)")
title("Abs of difference between real and predicted")
subtitle('pol = 1; test = 5; learning = 0.2; hidden = 256x(1,16); error = 1e-5')
colormap('turbo');
clim([min(min(abs(difference))), max(max(abs(difference)))])
colorbar()

xlabel("Position i in array")
ylabel("Position j in array")
