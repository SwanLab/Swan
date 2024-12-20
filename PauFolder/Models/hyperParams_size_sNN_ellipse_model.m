%% Testing with IRIS DATASET
clc;
clear;
close all;
%addpath ..;

% Remove path before running
addpath('..')
addpath('../Models')
addpath('../Datasets')
addpath('../../VariAuto/Codes')

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.2;
hiddenLayers    = 40.* ones(1, 6);

dLayers_array = 2.^(3:7);
nLayers_array = 2:8;
nReps = 3;

L2_error = zeros(length(dLayers_array), length(nLayers_array));
std_dev = zeros(length(dLayers_array), length(nLayers_array));


%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Chomog_ellipse.csv';

for id = 1:length(dLayers_array)
for in = 1:length(nLayers_array)

current_errors = zeros(1, nReps);

for ir = 1:nReps

hiddenLayers =  dLayers_array(id).* ones(1, nLayers_array(in));

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;

s.networkParams.costType = 'L2';
s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = [1, 2];
s.yFeatures = [3, 4, 5, 6];
cHomogIdxs = [11, 12, 22, 33];

% Load data
data   = cHomogData(s);
s.data = data;

% Train the model
opt = OptimizationProblem(s);
opt.solve();

%% Plot surface

% Load dataset from specified path
filePath = fullfile('..', 'Datasets', s.fileName);
tempData = readmatrix(filePath);

% Preallocate and evaluate y_data vector
yData = cell2mat(arrayfun(@(i) opt.computeOutputValues(tempData(i, 1:2)), 1:size(tempData, 1), 'UniformOutput',false)');

%% Compute error

Ytest_error = opt.computeOutputValues(data.Xtest) - data.Ytest;
L2test_error = zeros(1, size(Ytest_error, 2));
for i = 1:size(Ytest_error, 2)
    L2test_error(i) = norm(Ytest_error(:, i));
end

%current_errors(ir) = L2test_error(i);
current_errors(ir) = norm(L2test_error);

end

L2_error(id, in) = mean(current_errors);
std_dev(id, in) = std(current_errors);

end
end

%% Plot results

save('archi_errors.mat', 'L2_error', 'std_dev', 'nLayers_array', 'dLayers_array')

close all
hfig = figure;
tiledlayout(1, 2);

nexttile
surfer = surf(nLayers_array, dLayers_array, L2_error);
xlabel('Number of layers')
ylabel('Size of each layer')
zlabel('L2 error')

nexttile
surf(nLayers_array, dLayers_array, std_dev)
xlabel('Number of layers')
ylabel('Size of each layer')
zlabel('Std deviation of error')

adjust_figure_properties(hfig, font_size, picturewidth, hw_ratio)