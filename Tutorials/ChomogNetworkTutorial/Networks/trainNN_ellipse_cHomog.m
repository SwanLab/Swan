clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.2;
hiddenLayers    = 128 .* ones(1, 6);

%% INITIALIZATION 

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs     = 100; % 1000 is the best option, but we use 100 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

%% Get Training Data Ready
% I want to take advantage of physics
% behind problem to simplify it: cholesky decomposition as matrix will be
% SPD-> half of y features needed

load("KcoarseTraining.mat") %loads Kc object with cells with 8x8 training matrices
for i=1:length(Kc)
    Kc{i} = Kc{i} + 1e-10*eye(size(Kc{i})); % regularization -> (extremely small negative smallest eigenvalue (- 10^-16)
    Kc{i} = round(Kc{i},12); % if we don't round matrices are not SPD for some reason
end

processor = TrainingDataOrganizer(0.1:0.001:0.5,Kc);
trainData = processor.generateTrainingMatrix(); % Generate training data

%% Save Training Data as CSV
% Number of unique components from the symmetric 8x8 matrix
nFeatures = size(trainData,2) - 1;   % number of Y features (excluding Radius)

% First column is the input variable
colNames = ["Radius"];

% Generate stiffness component names
for i = 1:nFeatures
    colNames(end+1) = sprintf('Kcomp_%02d', i);
end

T = array2table(trainData, 'VariableNames', colNames);
writetable(T, 'trainData.csv'); % Save to CSV

s.fileName = 'trainData.csv';
%load("trainData.csv")
% Select the model's features --- 
s.xFeatures = 1;
s.yFeatures = 2:(size(trainData,2));
%cHomogIdxs = [11, 12, 22, 33];

% Load data
%data   = cHomogData(s);
data   = Data(s);
s.data = data;

% Train the model
opt = OptimizationProblemNN(s);
opt.solve();
opt.plotCostFnc();

save('Kcoarse_predictor.mat', 'opt'); % Save the model

%% Get Network Results

Xin = [0.25];

Y = opt.computeOutputValues(Xin);
dY = opt.computeGradient(Xin);
n = 8;
L = zeros(n);
idx = tril(true(n));
L(idx) = Y(:);                        % fill lower triangle
d = diag(L);
d(d <= 0) = eps;                      
L(1:n+1:end) = d;
K_coarse_pred = L*L.';%revert cholesky decomposition                       

fprintf('Output of the network:\n')
disp(K_coarse_pred)

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