%% Testing with IRIS DATASET
clc;
clear;
close all;
%addpath ..;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;  
lambda          = 0.0;
learningRate    = 0.1;
%hiddenLayers    = [8,20,20,20,8,6];
hiddenLayers    = [6, 24, 48, 48, 48, 24, 6, 3];

%% INITIALIZATION 
% try different feature combination, pairs of features enable the possibility of plotting boundaries
s.fileName = 'Chomog_ellipse.csv';
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;

% Select the model's features
s.xFeatures = [1, 2];
s.yFeatures = 3;

% Load data
data = cHomogData(s);

% Load model parameters
s.networkParams.hiddenLayers    = hiddenLayers;
s.data                          = data;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;

% Train the model
opt = OptimizationProblem(s);
opt.solve();

%% Plot surface

% Load dataset from specified path
filePath = fullfile('..', 'Datasets', s.fileName);
tempData = readmatrix(filePath);

% Preallocate and evaluate y_data vector
yData = arrayfun(@(i) opt.eval(tempData(i, 1:2)), 1:size(tempData, 1))';

% Set up figure for plotting
figure;

% Determine grid size for reshaping data
gridSize = floor(sqrt(size(tempData, 1)));

% Reshape data into grid format
gridA = reshape(tempData(:, 1), [gridSize, gridSize]);
gridB = reshape(tempData(:, 2), [gridSize, gridSize]);
gridC = reshape(tempData(:, 3), [gridSize, gridSize]);
gridY = reshape(yData, [gridSize, gridSize]);

% Plot 'Ground Truth' surface
surf(gridA, gridB, gridC, 'FaceColor', 'blue', 'EdgeColor', 'none', 'DisplayName', 'Ground Truth');
hold on;

% Plot 'Predicted' surface with transparency
surf(gridA, gridB, gridY, 'FaceColor', [1, 0.5, 0], 'EdgeColor', 'none', 'DisplayName', 'Predicted');
alpha(0.7);

% Add legend and axis labels
legend('Location', 'best');
xlabel('Sx');
ylabel('Sy');

% Display title for the specific tensor component
title('Component 00 of Constitutive Tensor');
