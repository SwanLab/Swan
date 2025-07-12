clc;
clear;
close all;

% Handle paths
% addpath(genpath('Tutorials'))
% addpath(genpath('src'))

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.0005;
hiddenLayers    = 40 .* ones(1, 10);

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Chomog_superformula_Big.csv';

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

% Select the model's features - x vector: [a, n] on n1 = n2 = n3
s.xFeatures = [1, 4];
s.yFeatures = [7, 8, 9, 10];
cHomogIdxs = [11, 12, 22, 33];

% Load data
%data   = cHomogData(s);
data   = Data(s);
s.data = data;

% Train the model
opt = OptimizationProblemNN(s);
opt.solve();
opt.plotCostFnc();

% Save the model
save('Tutorials/ChomogNetworkTutorial/SuperformNetwork.mat', 'opt')

%% Get Network Results

Xin = [0.25, 0.25];

Y = opt.computeOutputValues(Xin);
dY = opt.computeGradient(Xin);

fprintf('Output of the network:\n')
disp(Y)

fprintf('Jacobian of the network output w.r.t its input:\n')
disp(dY)

%% Compute error

Ytest_error = opt.computeOutputValues(data.Xtest) - data.Ytest;
L2test_error = zeros(1, size(Ytest_error, 2));
for i = 1:size(Ytest_error, 2)
    L2test_error(i) = norm(Ytest_error(:, i));
end

disp(norm(L2test_error))