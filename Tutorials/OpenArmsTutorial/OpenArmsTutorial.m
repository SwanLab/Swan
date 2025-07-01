clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 20;
lambda          = 0.0;
learningRate    = 0.2;
hiddenLayers    = 128 .* ones(1, 6);

%% INITIALIZATION 
% Store dataset file name
s.fileName = 'Resultados2.csv';

% Load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs = 1000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = [1, 2, 3, 4, 5, 6, 7];
s.yFeatures = [8];
% Load data
%data   = cHomogData(s);
%data = JuliaData(s);
data   = Data(s);
s.data = data;

% Train the model
opt = OptimizationProblemNN(s);
opt.solve();
opt.plotCostFnc();

% Obtain test data
[Xtest, Ytest] = opt.getTestData();  % Usa el método para obtener los datos de prueba

% Load the trained neural network
network = opt.getNetwork();

%{
% Initialize Ypred as an empty column vector
Ypred = zeros(size(Ytest));

% Pass every Xtest sample to the neural network and store predictions
for i = 1:size(Xtest, 1)
    %Ypred(i) = network.forwardprop(Xtest(i, :), Ytest(i, :));  % Prediction for each case
    Ypred(i) = network.computeYOut(Xtest(i, :));
end
%}
Ypred = network.computeYOut(Xtest); % Vectorize the computation of Ypred

%(debugging) Histogram for the distribution of Ypred
figure;
edges = linspace(-1, 2, 31);  % 30 bins between -1 and 2
histogram(Ypred, edges);
%histogram(Ypred, 'NumBins', 10); 
title('Distribution of predicted Ytest');

figure;
edges = linspace(-1, 2, 31);  % 30 bins between -1 and 2
histogram(Ytest, edges);
%histogram(Ytest, 'NumBins', 10);
title('Distribution of Test Y');

% Denormalization
Xtest = Xtest .* data.sigmaX + data.muX;
Ypred = Ypred .* data.sigmaY + data.muY;
Ytest = Ytest .* data.sigmaY + data.muY;

% Consumption dependance on speed cubed
figure;
plot(Xtest(:,4),Ytest,'o')
xlabel('Speed cubed (m/s)^3')
ylabel('Fuel consumption')


% Compute the mean square error (MSE)
mse = mean((Ypred - Ytest).^2);

disp(['Error cuadrático medio (MSE) en los datos de prueba: ', num2str(mse)]);

% Compute the difference between real and predicted values
difference = Ytest - Ypred;

% Convert Xtest to a table and rename the columns
input_data = array2table(Xtest);
input_data.Properties.VariableNames =  {'rpm','Windy(cosine)', 'Windy(m/s)', 'Speed^3 (m/s)^3', 'Yaw', 'Pitch', 'Roll'};  

% Create a table of outputs and rename the columns
output_data = table(Ytest, Ypred, difference);
output_data.Properties.VariableNames = {'Cons. real', 'Cons. prediction', 'Difference'};

% Combine input and output tables
result_table = [input_data, output_data];

% Display the resulting table
disp('Tabla de resultados:');
disp(result_table);
