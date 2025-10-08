clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;
lambda          = 0.0;
learningRate    = 0.001;
hiddenLayers    = 256 .* ones(1, 8);  % Capas más profundas para matriz 8x8

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'DB_stiffness_matrix.csv';  % Nuevo dataset

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs     = 1000; % Más épocas para matriz compleja
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
% X: 1 feature (parámetro de entrada)
% Y: 64 features (matriz 8x8 aplanada)
s.xFeatures = [1];  % Una sola característica de entrada
s.yFeatures = 2:65; % 64 características de salida (8x8 = 64)

% Load data
data   = Data(s);
s.data = data;

% Train the model
opt = OptimizationProblemNN(s);
opt.solve();
opt.plotCostFnc();

% Save the model
save('Tutorials/ChomogNetworkTutorial/Networks/network_stiffness_matrix.mat', 'opt')

%% Get Network Results

Xin = [0.5];  % Valor de entrada de ejemplo

Y = opt.computeOutputValues(Xin);
dY = opt.computeGradient(Xin);

% Reshape output to 8x8 matrix
stiffness_matrix = reshape(Y, 8, 8);

fprintf('Output of the network (8x8 stiffness matrix):\n')
disp(stiffness_matrix)

fprintf('Jacobian of the network output w.r.t its input:\n')
disp(dY)

%% Plot results

% Load dataset from specified path
filePath = fullfile('Tutorials', 'ChomogNetworkTutorial', 'Datasets', s.fileName);
tempData = readmatrix(filePath);

% Preallocate and evaluate y_data vector
yData = cell2mat(arrayfun(@(i) opt.computeOutputValues(tempData(i, 1)), 1:size(tempData, 1), 'UniformOutput',false)');

% Plot some components of the stiffness matrix
figure;
tiledlayout(2, 2)

% Plot diagonal elements
for i = 1:4
    nexttile
    scatter(tempData(:, 1), tempData(:, s.yFeatures(i)), 'b', 'filled', 'DisplayName', 'Ground Truth');
    hold on;
    scatter(tempData(:, 1), yData(:, i), 'r', 'filled', 'DisplayName', 'Predicted');
    xlabel('Input Parameter');
    ylabel(['Stiffness Matrix Element (', num2str(i), ',', num2str(i), ')']);
    legend('Location', 'best');
    title(['Diagonal Element K(', num2str(i), ',', num2str(i), ')']);
end

% Adjust figure size and position
gcf.Position = [100 100 1000 600];

%% Compute error

Ytest_error = opt.computeOutputValues(data.Xtest) - data.Ytest;
L2test_error = zeros(1, size(Ytest_error, 2));
for i = 1:size(Ytest_error, 2)
    L2test_error(i) = norm(Ytest_error(:, i));
end

fprintf('L2 Test Error: %f\n', norm(L2test_error))

%% Test matrix properties
% Verificar que la matriz de rigidez sea simétrica y definida positiva
fprintf('\nMatrix Properties:\n')
fprintf('Is symmetric: %s\n', mat2str(issymmetric(stiffness_matrix, 1e-10)))
fprintf('Is positive definite: %s\n', mat2str(all(eig(stiffness_matrix) > 0)))
fprintf('Condition number: %f\n', cond(stiffness_matrix))
