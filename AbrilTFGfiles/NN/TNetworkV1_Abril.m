% This script is intended to train the NN for the Downscaling dataset V1

clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 20;
lambda          = 0.0;
learningRate    = 0.05;
hiddenLayers    = 6 .* ones(1,2);
 

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Tdata.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 1000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the T column to train
T_type='T3'

% Select the model's features
s.xFeatures = 1:3;

switch T_type
    case 'T1'
        s.yFeatures = [4,5];      %T1
    case 'T2'         
        s.yFeatures = [6,7];      %T2
    case 'T3'         
        s.yFeatures = [8,9];      %T3
    case 'T4'   
        s.yFeatures = [10,11];    %T4
    case 'T5'   
        s.yFeatures = [12,13];    %T5
    case 'T6'   
        s.yFeatures = [14,15];    %T6
    case 'T7'
        s.yFeatures = [16,17];    %T7
    case 'T8'   
        s.yFeatures = [18,19];    %T8
end


% Load data
%data   = cHomogData(s);
data   = Data(s);
s.data = data;

% Train the model
T_NN = OptimizationProblemNN(s);
T_NN.solve();
T_NN.plotCostFnc();
MSETrain    = immse(T_NN.computeOutputValues(data.Xtrain), data.Ytrain);


%% Plot surface

% Load dataset from specified path
filePath = fullfile('AbrilTFGfiles', s.fileName);
tempData = readmatrix(filePath);

real = tempData(:, s.yFeatures);
predicted = zeros(size(real));

for i = 1:size(real,1)
    predicted(i, :) = T_NN.computeOutputValues(tempData(i, s.xFeatures));
end

difference = real-predicted;
disp(max(abs(difference)));

switch T_type
    case 'T1'
        string ="T1.mat"; %T1
    case 'T2'
        string ="T2.mat"; %T2
    case 'T3'
        string ="T3.mat"; %T3
    case 'T4'
        string ="T4.mat"; %T4
    case 'T5'
        string ="T5.mat"; %T5
    case 'T6'
        string ="T6.mat"; %T6
    case 'T7'
        string ="T7.mat"; %T7
    case 'T8'
        string ="T8.mat"; %T8
end


FileName=fullfile('AbrilTFGfiles','NN',string);
    save(FileName, "T_NN","real","predicted","difference","T_type");


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
