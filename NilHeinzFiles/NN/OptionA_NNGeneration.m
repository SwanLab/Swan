%This script builds and trains 3 neural networks:
% - KcoarsePredictorNN : For a radius input, predicts Coarse Stiffness
% matrix Kcoasrse
% - McoarsePredictorNN : For a radius input, predicts Coarse Mass
% matrix Mcoarse
% - TPredictorNN : FOr r,(x,y) input, predicts T value. Used for Method A
% to regenerate T for different subdomains.

close all;

% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 20;
lambda          = 0.0;
learningRate    = 0.2;
hiddenLayers    = 128 .* ones(1, 6);

%% Kcoarse Predict
s1.polynomialOrder = pol_deg;
s1.testRatio       = testratio;
s1.networkParams.hiddenLayers    = hiddenLayers;
s1.optimizerParams.learningRate  = learningRate;
s1.optimizerParams.maxEpochs     = 5000; % 1000 is the best option, but we use 100 to pass the tutorial quickly
s1.costParams.lambda             = lambda;
s1.costParams.costType           = 'L2';

s1.networkParams.HUtype = 'ReLU';
s1.networkParams.OUtype = 'linear';

s1.fileName = 'KcoarseTrainingData.csv';
s1.xFeatures = 1;
s1.yFeatures = 2:37;

% Load Data
data   = Data(s1);
s1.data = data;

% Train the model
KcPredictorNN = OptimizationProblemNN(s1);
KcPredictorNN.solve();
KcPredictorNN.plotCostFnc();

save('KcoarsePredictorNN.mat', 'KcPredictorNN'); % Save the model. 


%% Mcoarse Predict
s2.polynomialOrder = pol_deg;
s2.testRatio       = testratio;
s2.networkParams.hiddenLayers    = hiddenLayers;
s2.optimizerParams.learningRate  = learningRate;
s2.optimizerParams.maxEpochs     = 2000; % 1000 is the best option, but we use 100 to pass the tutorial quickly
s2.costParams.lambda             = lambda;
s2.costParams.costType           = 'L2';

s2.networkParams.HUtype = 'ReLU';
s2.networkParams.OUtype = 'linear';

s2.fileName = 'McTrainingData.csv';
s2.xFeatures = 1;
s2.yFeatures = 2:37;

% Load Data
data   = Data(s2);
s2.data = data;

% Train the model
McPredictorNN = OptimizationProblemNN(s2);
McPredictorNN.solve();
McPredictorNN.plotCostFnc();

save('McoarsePredictorNN.mat', 'McPredictorNN'); % Save the model


%% Method A- Direct NN to predict full T

% INITIALIZATION - load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs     = 500; 
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

s.fileName = 'TTrainingData.csv';
s.xFeatures = 1:3;
s.yFeatures = 4:19;

% Load Data
data   = Data(s);
s.data = data;

% Train the model
tPredictorNN = OptimizationProblemNN(s);
tPredictorNN.solve();
tPredictorNN.plotCostFnc();

save('TpredictorNN.mat', 'tPredictorNN'); % Save the model
Xin = [0.25, 0.1,0.2];

Y = tPredictorNN.computeOutputValues(Xin);

