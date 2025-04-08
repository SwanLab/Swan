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
lambda_array    = linspace(0, 0.01, 10);
learningRate    = 0.2;
hiddenLayers    = 12.* ones(1, 3);

nReps = 1;

L2_test_error = zeros(size(lambda_array));
std_test_dev =  zeros(size(lambda_array));

L2_train_error =  zeros(size(lambda_array));
std_train_dev =  zeros(size(lambda_array));


%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Chomog_ellipse.csv';

for il = 1:size(lambda_array, 2)

current_test_errors = zeros(1, nReps);
current_train_errors = zeros(1, nReps);

for ir = 1:nReps

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda_array(il);

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
current_test_error = zeros(1, size(Ytest_error, 2));
for i = 1:size(Ytest_error, 2)
    current_test_error(i) = norm(Ytest_error(:, i));
end

Ytrain_error = opt.computeOutputValues(data.Xtrain) - data.Ytrain;
current_train_error = zeros(1, size(Ytrain_error, 2));
for i = 1:size(Ytrain_error, 2)
    current_train_error(i) = norm(Ytrain_error(:, i));
end

current_test_errors(ir) = norm(current_test_error);
current_train_errors(ir) = norm(current_train_error);



end

L2_test_error(il) = mean(current_test_errors);
std_test_dev(il) = std(current_test_errors);

L2_train_error(il) = mean(current_train_errors);
std_train_dev(il) = std(current_train_errors);

end

%% Plot results

save('lambda_errors.mat', 'L2_test_error', 'std_test_dev', 'L2_train_error', 'std_train_dev', 'lambda_array')

close all
hfig = figure;
%tiledlayout(1, 2);

%nexttile
grid on
hold on
plot(lambda_array, L2_test_error, '-b', 'Linewidth', 1.25, 'DisplayName', 'Test error');
plot(lambda_array, L2_train_error, '-r', 'Linewidth', 1.25, 'DisplayName', 'Training error');
hold off
xlabel('$\lambda$', 'Interpreter','latex')
ylabel('L2 error', 'Interpreter','latex')
legend('Location', 'best')

%nexttile
%hold on
%plot(lambda_array, std_test_dev);
%plot(lambda_array, std_train_dev);
%hold off
%xlabel('$\lambda$', 'Interpreter','latex')
%ylabel('Std deviation error', 'Interpreter','latex')

adjust_figure_properties(hfig, 16, 25, 0.47)