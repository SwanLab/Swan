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
hiddenLayers    = 128.* ones(1, 5);
inLayer = 2;
outLayer = 4;
nTaper = 1:2;
altLayers(1,:) = hiddenLayers;
altLayers(2,:) = [floor((hiddenLayers(1)+inLayer) / 2), hiddenLayers(2:end-1), floor((hiddenLayers(1)+outLayer) / 2)];
altLayers(3,:) = [floor(inLayer+(hiddenLayers(1)-inLayer)/3), floor(inLayer+2*(hiddenLayers(1)-inLayer)/3), hiddenLayers(3:end-2),floor(outLayer+2*(hiddenLayers(1)-inLayer)/3),floor(outLayer+(hiddenLayers(1)-inLayer)/3)];

dLayers_array = 2.^(3:7);
nLayers_array = 2:8;
nReps = 24;

L2_error = zeros(size(altLayers, 1), 1);
std_dev = zeros(size(altLayers, 1), 1);


%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Chomog_ellipse.csv';

for it = 1:size(altLayers, 1)

current_errors = zeros(1, nReps);

for ir = 1:nReps

hiddenLayers =  altLayers(it, :);

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

L2_error(it) = mean(current_errors);
std_dev(it) = std(current_errors);

end

%% Plot results

save('taper_errors.mat', 'L2_error', 'std_dev', 'altLayers')

close all
hfig = figure;
tiledlayout(1, 2);

nexttile
bar(0:size(altLayers, 1)-1, L2_error);
xlabel('Taper case', 'Interpreter','latex')
ylabel('L2 error', 'Interpreter','latex')

nexttile
bar(0:size(altLayers, 1)-1, std_dev);
xlabel('Taper case', 'Interpreter','latex')
ylabel('Std deviation error', 'Interpreter','latex')

%adjust_figure_properties(hfig, font_size, picturewidth, hw_ratio)