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


%% Función para hacer predicciones en nuevas situaciones
function predicciones = predecirValor(inputData)
    % inputData es un vector de situaciones que queremos predecir
    % Asegúrate de que la estructura 'opt' y el modelo de regresión
    % están disponibles, en caso de que desees usar el modelo entrenado.
    
    % Predecir usando la red entrenada
    predicciones = opt.predict(inputData); 
end

%% Ejemplo de uso: predicción para nuevas situaciones
nuevasSituaciones = [5, 10;   % Ejemplo de entrada: [x1, x2]
                     2, 8];   % Ejemplo de entrada: [x1, x2]

% Llamada a la función de predicción
resultados = predecirValor(nuevasSituaciones);

% Mostrar los resultados
disp('Predicciones para las nuevas situaciones:');
disp(resultados);