clc
clear
close all

rng(1)

fprintf('\n')
fprintf('============================================\n')
fprintf('      LATTICE NN SENSITIVITY LOSS CHECK\n')
fprintf('============================================\n')


%% ------------------------------------------------------------
%  DATA
% -------------------------------------------------------------

sData.fileName = ...
    'LatticeNNData.csv';

sData.polynomialOrder = ...
    3;

sData.testRatio = ...
    30;

sData.xFeatures = ...
    1:3;

sData.yFeatures = ...
    4:9;

data = ...
    Data(sData);


fprintf('\nMain C database loaded\n')

fprintf(...
    'Training samples = %d\n',...
    size(data.Xtrain,1))

fprintf(...
    'Test samples     = %d\n',...
    size(data.Xtest,1))


%% ------------------------------------------------------------
%  NETWORK
%
%  Same architecture as the current baseline.
% -------------------------------------------------------------

s.networkParams.hiddenLayers = ...
    [100 200 300 200 100 50];

s.networkParams.HUtype = ...
    'SoftPlus';

s.networkParams.OUtype = ...
    'linear';


%% ------------------------------------------------------------
%  COST
% -------------------------------------------------------------

s.costParams.costType = ...
    'L2';

s.costParams.lambda = ...
    1e-9;


%% ------------------------------------------------------------
%  OPTIMIZER
%
%  We are NOT training in this tutorial.
%  These values are only supplied because the optimization
%  problem creates its optimizer during construction.
% -------------------------------------------------------------

s.optimizerParams.learningRate = ...
    5e-3;

s.optimizerParams.maxEpochs = ...
    1;

s.optimizerParams.earlyStop = ...
    1;


%% ------------------------------------------------------------
%  SENSITIVITY DATABASE
% -------------------------------------------------------------

s.sensitivityParams.sensitivityFileName = ...
    'LatticeNNSensitivityD4.csv';


%% ------------------------------------------------------------
%  dC / drho
% -------------------------------------------------------------

s.sensitivityParams.variables(1).name = ...
    'rho';

s.sensitivityParams.variables(1).featureIndex = ...
    1;

s.sensitivityParams.variables(1).targetColumns = ...
    4:9;

s.sensitivityParams.variables(1).h = ...
    1e-5;

s.sensitivityParams.variables(1).weight = ...
    1;


%% ------------------------------------------------------------
%  dC / da
% -------------------------------------------------------------

s.sensitivityParams.variables(2).name = ...
    'a';

s.sensitivityParams.variables(2).featureIndex = ...
    2;

s.sensitivityParams.variables(2).targetColumns = ...
    10:15;

s.sensitivityParams.variables(2).h = ...
    9.8e-6;

s.sensitivityParams.variables(2).weight = ...
    1;


%% ------------------------------------------------------------
%  dC / db
% -------------------------------------------------------------

s.sensitivityParams.variables(3).name = ...
    'b';

s.sensitivityParams.variables(3).featureIndex = ...
    3;

s.sensitivityParams.variables(3).targetColumns = ...
    16:21;

s.sensitivityParams.variables(3).h = ...
    2.7e-5;

s.sensitivityParams.variables(3).weight = ...
    1;


%% ------------------------------------------------------------
%  COMPLETE PROBLEM DATA
% -------------------------------------------------------------

s.data = ...
    data;


%% ------------------------------------------------------------
%  CREATE PROBLEM
%
%  IMPORTANT:
%  ReLoBRaLo is NOT enabled yet.
% -------------------------------------------------------------

problemNN = ...
    OptimizationProblemNN(s);


fprintf('\nOptimization problem created.\n')
fprintf('Three sensitivity functionals enabled.\n')


%% ------------------------------------------------------------
%  SENSITIVITY DIAGNOSTICS
% -------------------------------------------------------------

problemNN.printSensitivityDiagnostics();


fprintf('\n')
fprintf('============================================\n')
fprintf('      SENSITIVITY LOSS CHECK COMPLETED\n')
fprintf('============================================\n')