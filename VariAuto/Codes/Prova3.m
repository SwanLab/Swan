%% Testing with Regression on Custom Data
clc;
clear;
close all;

%% Inicialización de hiperparámetros
pol_deg         = 1;                
testratio       = 20;              
lambda          = 0;                
learningRate    = 0.1;              
hiddenLayers    = [3,2,1,2,3];      

%% INICIALIZACIÓN
s.fileName = '/OpenArms/archivo_modificado.csv';          
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.xFeatures = [1, 2, 3, 4, 5, 6, 7];                
s.yFeatures = [8];                 
data = Data(s);

s.networkParams.costType = 'L2';
s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';
s.networkParams.hiddenLayers    = hiddenLayers;
s.data                          = data;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;

opt = OptimizationProblem(s);
opt.solve();

% Obtener los datos de prueba
[Xtest, Ytest] = opt.getTestData();  % Usa el método para obtener los datos de prueba

% Realizar predicciones con la red neuronal entrenada
network = opt.getNetwork();

% Inicializamos Ypred como un columna vacía
Ypred = zeros(size(Ytest));

% Pasamos cada muestra de Xtest a la red neuronal y guardamos las predicciones
for i = 1:size(Xtest, 1)
    Ypred(i) = network.forwardprop(Xtest(i, :), Ytest(i, :));  % Predicción para cada caso
end

% Calcular el error cuadrático medio (MSE)
mse = mean((Ypred - Ytest).^2);  
disp(['Error cuadrático medio (MSE) en los datos de prueba: ', num2str(mse)]);

% Crear la diferencia entre los valores reales y predichos
difference = Ytest - Ypred;

% Convertir Xtest en una tabla y renombrar las columnas
input_data = array2table(Xtest);  % Convertir Xtest en tabla
input_data.Properties.VariableNames = {'Rpm', 'Windy(º)', 'Windy(m/s)', 'Speed (m/s)', 'Yaw', 'Pitch', 'Roll'};  % Renombrar columnas

% Crear la tabla de salidas y renombrar las columnas
output_data = table(Ytest, Ypred, difference);  % Crear la tabla de salidas
output_data.Properties.VariableNames = {'Cons. real', 'Cons. prediction', 'Difference'};  % Renombrar columnas

% Combinar las tablas de entradas y salidas
result_table = [input_data, output_data];

% Mostrar la tabla final
disp('Tabla de resultados:');
disp(result_table);



