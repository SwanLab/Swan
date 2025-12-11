clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 8;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.01;
hiddenLayers    =[16 20 16 6 3];

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'DataQ.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 20000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = [1];
cHomogIdxs = [11, 12, 22, 33];

% Load data
%data   = cHomogData(s);


% Train the model

for i=1:15
    s.yFeatures = [i+1];
    data   = Data(s);
    s.data = data;
    Q_NN{i} = OptimizationProblemNN(s);
    Q_NN{i}.solve();
    Q_NN{i}.plotCostFnc();

    %MSETrain    = immse(Q_NN(i).computeOutputValues(data.Xtrain), data.Ytrain);
end

string ="Q_NN.mat";
FileName=fullfile('AbrilTFGfiles','NN',string);
    save(FileName, "Q_NN");

%% Plot surface

% Load dataset from specified path
filePath = fullfile('AbrilTFGfiles', s.fileName);
tempData = readmatrix(filePath);

real = tempData(:,s.yFeatures);
predicted = zeros(size(real));

for i = 1:size(real,1)
    predicted(i,:) = Q_NN.computeOutputValues(tempData(i,s.xFeatures));
end

difference = real-predicted;
disp(norm(difference));
