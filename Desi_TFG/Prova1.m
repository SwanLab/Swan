%% Testing with Regression on Custom Data
clc;
clear;
close all;


%% Initialization of hyperparameters
pol_deg         = 1;                
testratio       = 30;               
lambda          = 0;                
learningRate    = 0.1;              
hiddenLayers    = [3,2,1,2,3];      

%% INITIALIZATION
s.fileName = '/OpenArms/Resultados2.csv';          
s.polynomialOrder = pol_deg;        %% Orden del polinomio
s.testRatio       = testratio;      %% Proporción de datos de prueba
s.xFeatures = [1, 2];
s.yFeatures = [1];
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
%opt.plotRegressionResults();


% if data.nFeatures == 2  % If you want to be asked for Features change it in "Data" Class
%     opt.plotRegressionBoundary();
% end

% Datos para evaluar
x_new = [1600, 133, 2, 2.5, 0.75, 1.25, -1];  % 7 características

% x_new a matriz fila
x_new = reshape(x_new, 1, []);

% Predición
y_pred = opt.computeOutputValues(x_new);  % Usamos el método para obtener la salida
disp(['Predicción: ', num2str(y_pred)]);




