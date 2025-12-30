clc;
clear;
close all;

%% Case parameters
p.nelem      = 20;
p.Sampling   ='Isolated';     %'Isolated'/'Oversampling'
p.Inclusion  ='HoleRaul';    %'Material'/'Hole'/'HoleRaul
meshName     = p.nelem+"x"+p.nelem;
%% Initialization of hyperparameters
pol_deg         = 6;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.001;
hiddenLayers    =[100 150 200 250 200 150 100 60 40 24 20];

%% INITIALIZATION 
% Store dataset file name
s.fileName = fullfile('AbrilTFGfiles',"Data",p.Inclusion,p.Sampling,meshName,'DataQ.csv');

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 10000000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = [1];
s.yFeatures = [2:1:11];
cHomogIdxs = [11, 12, 22, 33];

% Load data
data   = Data(s);
s.data = data;

% Train the model

Q_NN = OptimizationProblemNN(s);
Q_NN.solve();
Q_NN.plotCostFnc();


string ="Q_NN.mat";
FileName=fullfile('AbrilTFGfiles',"Data",p.Inclusion,p.Sampling,string);
save(FileName, "Q_NN","pol_deg");

%% Plot surface

% Load dataset from specified path
filePath = s.fileName;
tempData = readmatrix(filePath);

real = tempData(:,s.yFeatures);
predicted = zeros(size(real));

%for i = 1:size(real,1)
%    predicted(i,:) = Q_NN.computeOutputValues(tempData(i,s.xFeatures));
%end
%
%difference = real-predicted;
%disp(norm(difference));
