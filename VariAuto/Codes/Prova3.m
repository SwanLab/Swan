%% Testing with Regression on Custom Data
clc;
clear;
close all;
format compact

%% Inicialización de hiperparámetros
pol_deg         = 1;                
testratio       = 20;              
lambda          = 0.01;                
learningRate    = 0.2;              
hiddenLayers    = [10,5,10];  %[3,2,1,2,3];
% tanh

%% INICIALIZACIÓN
s.fileName = 'OpenArms/Resultados2.csv'; 
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.xFeatures = [1, 2, 3, 4, 5, 6, 7];  
s.nLayers = length(hiddenLayers);
s.neuronsPerLayer = hiddenLayers;
s.yFeatures = [8];                 
data = Data(s);
s.networkParams.costType = 'L2';
s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';
%s.networkParams.HUtype = 'tanh';
%s.networkParams.OUtype = 'linear';
s.networkParams.hiddenLayers    = hiddenLayers;
s.data                          = data;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;
s.costParams.costType = 'L2';

opt = OptimizationProblem(s);
learnableVariables = LearnableVariables(s);
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


%{
%(debugging) Histograma para ver la distribución de Ypred, [0,1]
figure;
edges = linspace(-0.05, 1.05, 23);  % 10 bins between 0 and 1
histogram(Ypred, edges);
%histogram(Ypred, 'NumBins', 10); 
title('Distribution of predicted Ytest');

figure;
edges = linspace(0, 1, 21);  % 10 bins between 0 and 1
histogram(Ytest, edges);
%histogram(Ytest, 'NumBins', 10);
title('Distribution of Test Y');
%}

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

%{
Xtest = Xtest .* (data.Xmax - data.Xmin) + data.Xmin;
Ypred = Ypred .* (data.Ymax - data.Ymin) + data.Ymin;
Ytest = Ytest .* (data.Ymax - data.Ymin) + data.Ymin;
%}

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



