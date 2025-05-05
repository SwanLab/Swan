clc;
clear;
close all;

% Handle paths
addpath('src/NeuralNetwork')
addpath('src/Problems/Optimization')
addpath('Tutorials/ChomogNetworkTutorial')
addpath('Tutorials/ChomogNetworkTutorial/Datasets')

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;
lambda          = 0.000001;
learningRate    = 0.1;
hiddenLayers    = 6 .* ones(1, 2);

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Plane_DB.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = [1, 2];
s.yFeatures = [3];

% Jacobian solution
xSlope = 0.2;
ySlope = 3.9;

% Load data
%data   = cHomogData(s);
data   = Data(s);
s.data = data;

% Train the model
opt = OptimizationProblem(s);
opt.solve();
opt.plotCostFnc();

%% Get Network Results

Xin = [0.25, 0.25];

Y = opt.computeOutputValues(Xin);
dY = opt.computeGradient(Xin);

fprintf('Output of the network:\n')
disp(Y)

fprintf('Jacobian of the network output w.r.t its input:\n')
disp(dY)

%% Plot surface

% Load dataset from specified path
filePath = fullfile('Tutorials', 'ChomogNetworkTutorial', 'Datasets', s.fileName);
tempData = readmatrix(filePath);

% Preallocate and evaluate y_data vector
yData = cell2mat(arrayfun(@(i) opt.computeOutputValues(tempData(i, 1:2)), 1:size(tempData, 1), 'UniformOutput',false)');
dydaData = cell2mat(arrayfun(@(i) dot([1, 0], opt.computeGradient(tempData(i, 1:2))), 1:size(tempData, 1), 'UniformOutput',false)');
dydbData = cell2mat(arrayfun(@(i) dot([0, 1], opt.computeGradient(tempData(i, 1:2))), 1:size(tempData, 1), 'UniformOutput',false)');

% Determine grid size for reshaping data
gridSize = floor(sqrt(size(tempData, 1)));

% Reshape data into grid format
gridA = reshape(tempData(:, 1), [gridSize, gridSize]);
gridB = reshape(tempData(:, 2), [gridSize, gridSize]);

% Set up plot
hfig = figure;
tiledlayout(1, 3)

% Plot output (predicted vs real)
nexttile
% Grid the output values
gridC = reshape(tempData(:, s.yFeatures(1)), [gridSize, gridSize]);
gridY = reshape(yData(:, 1), [gridSize, gridSize]);

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
title('Network output');

% Plot Jacobian (predicted vs real)
nexttile
% Grid the output values
gridC = reshape(tempData(:, s.yFeatures(1)).*0 + ySlope, [gridSize, gridSize]);
gridY = reshape(dydbData(:, 1), [gridSize, gridSize]);

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
title('X slope');


% Plot Jacobian (predicted vs real)
nexttile
% Grid the output values
gridC = reshape(tempData(:, s.yFeatures(1)).*0 + xSlope, [gridSize, gridSize]);
gridY = reshape(dydaData(:, 1), [gridSize, gridSize]);

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
title('Y slope');

% Adjust figure size and position
hfig.Position = [100 100 1000 600];

%% Compute error

Ytest_error = opt.computeOutputValues(data.Xtest) - data.Ytest;
L2test_error = zeros(1, size(Ytest_error, 2));
for i = 1:size(Ytest_error, 2)
    L2test_error(i) = norm(Ytest_error(:, i));
end

disp(norm(L2test_error))