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
s.optimizerParams.maxEpochs = 3000; % 1000 is the best option, but we use 10 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

% Select the model's features
s.xFeatures = [1, 2, 3, 4, 5, 6, 7];
s.yFeatures = [8];

% Load data
%data   = cHomogData(s);
data   = Data(s);
s.data = data;

% Train the model
opt = OptimizationProblemNN(s);
opt.solve();
opt.plotCostFnc();

% Obtener los datos de prueba
[Xtest, Ytest] = opt.getTestData();  % Usa el método para obtener los datos de prueba

% Realizar predicciones con la red neuronal entrenada
network = opt.getNetwork();

% Inicializamos Ypred como un columna vacía
Ypred = zeros(size(Ytest));

% Pasamos cada muestra de Xtest a la red neuronal y guardamos las predicciones
for i = 1:size(Xtest, 1)
    %Ypred(i) = network.forwardprop(Xtest(i, :), Ytest(i, :));  % Predicción para cada caso
    Ypred(i) = network.computeYOut(Xtest(i, :));
end

%(debugging) Histograma para ver la distribución de Ypred, [-1,1]
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


% Calcular el error cuadrático medio (MSE)
mse = mean((Ypred - Ytest).^2);

disp(['Error cuadrático medio (MSE) en los datos de prueba: ', num2str(mse)]);

% Crear la diferencia entre los valores reales y predichos
difference = Ytest - Ypred;


% Convertir Xtest en una tabla y renombrar las columnas
input_data = array2table(Xtest);  % Convertir Xtest en tabla
input_data.Properties.VariableNames =  {'rpm','Windy(cosine)', 'Windy(m/s)', 'Speed^3 (m/s)^3', 'Yaw', 'Pitch', 'Roll'};  %'rpm', Renombrar columnas

% Crear la tabla de salidas y renombrar las columnas
output_data = table(Ytest, Ypred, difference);  % Crear la tabla de salidas
output_data.Properties.VariableNames = {'Cons. real', 'Cons. prediction', 'Difference'};  % Renombrar columnas

% Combinar las tablas de entradas y salidas
result_table = [input_data, output_data];

% Mostrar la tabla final
disp('Tabla de resultados:');
disp(result_table);

%% Plot surface
%{
% Load dataset from specified path
filePath = fullfile('Tutorials', 'ChomogNetworkTutorial', s.fileName);
tempData = readmatrix(filePath);

% Preallocate and evaluate y_data vector
yData = cell2mat(arrayfun(@(i) opt.computeOutputValues(tempData(i, 1:2)), 1:size(tempData, 1), 'UniformOutput',false)');

Y = opt.computeOutputValues([0.25, 0.25]);
dY = opt.computeGradient([0.25, 0.25]);

disp(Y)
disp(dY)

Y = opt.computeOutputValues([0.27, 0.27]);
dY = opt.computeGradient([0.27, 0.27]);

disp(Y)
disp(dY)

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
%}