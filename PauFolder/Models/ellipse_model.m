%% Testing with IRIS DATASET
clc;
clear;
close all;
%addpath ..;

% Remove path before running
%rmpath(genpath(pwd))
addpath('..')
addpath('../Models')
addpath('../Datasets')
addpath('../../VariAuto/Codes')

%% Initialization of hyperparameters
pol_deg         = 3;
testratio       = 30;  
lambda          = 0.00001;
learningRate    = 0.15;
%hiddenLayers    = [8,10,10,10,8,6];
%hiddenLayers    = [6, 24, 48, 48, 48, 24, 6, 3];
hiddenLayers = [32, 32, 16];

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Chomog_ellipse.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;

% Select the model's features
s.xFeatures = [1, 2];
yFeaturesArray = 3:6;
cHomogIdxs = [11, 12, 22, 33];

% Set up plot
hfig = figure;
tiledlayout(2, 2)

for i = 1:length(yFeaturesArray)

    % Select the model's Y features
    s.yFeatures = yFeaturesArray(i);
    
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
    yData = arrayfun(@(i) opt.eval(tempData(i, 1:2)), 1:size(tempData, 1))';
    
    % Set up figure for plotting
    nexttile
    %figure;
    
    % Determine grid size for reshaping data
    gridSize = floor(sqrt(size(tempData, 1)));
    
    % Reshape data into grid format
    gridA = reshape(tempData(:, 1), [gridSize, gridSize]);
    gridB = reshape(tempData(:, 2), [gridSize, gridSize]);
    gridC = reshape(tempData(:, s.yFeatures), [gridSize, gridSize]);
    gridY = reshape(yData, [gridSize, gridSize]);
    
    % Plot 'Ground Truth' surface
    surf(gridA, gridB, gridC, 'FaceColor', 'blue', 'EdgeColor', 'none', 'DisplayName', 'Ground Truth');
    hold on;
    
    % Plot 'Predicted' surface with transparency
    surf(gridA, gridB, gridY, 'FaceColor', [1, 0.5, 0], 'EdgeColor', 'none', 'DisplayName', 'Predicted');
    alpha(0.7);
    view(125,20)
    
    % Add legend and axis labels
    legend('Location', 'best');
    xlabel('Sx');
    ylabel('Sy');
    
    % Display title for the specific tensor component
    title(['Component ', num2str(cHomogIdxs(i)),' of Constitutive Tensor']);

end
% Adjust figure size and position
hfig.Position = [100 100 1000 600];
