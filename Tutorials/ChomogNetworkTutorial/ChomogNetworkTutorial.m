clc;
clear;
close all;

% Handle paths
addpath(genpath('Tutorials'))
addpath(genpath('src'))

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.2;
hiddenLayers    = 128 .* ones(1, 6);

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Chomog_ellipse.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 100; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = [1, 2];
s.yFeatures = [3, 4, 5, 6];
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
save('Tutorials/ChomogNetworkTutorial/ChomogNetwork.mat', 'opt')

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

% Determine grid size for reshaping data
gridSize = floor(sqrt(size(tempData, 1)));

% Reshape data into grid format
gridA = reshape(tempData(:, 1), [gridSize, gridSize]);
gridB = reshape(tempData(:, 2), [gridSize, gridSize]);

% Set up plot
hfig = figure;
tiledlayout(2, 2)

for i = 1:length(s.yFeatures)

    gridC = reshape(tempData(:, s.yFeatures(i)), [gridSize, gridSize]);
    gridY = reshape(yData(:, i), [gridSize, gridSize]);
    
    % Plot 'Ground Truth' surface
    nexttile
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

%% Compute error

Ytest_error = opt.computeOutputValues(data.Xtest) - data.Ytest;
L2test_error = zeros(1, size(Ytest_error, 2));
for i = 1:size(Ytest_error, 2)
    L2test_error(i) = norm(Ytest_error(:, i));
end

disp(norm(L2test_error))