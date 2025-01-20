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
s.fileName = 'Iris.csv';          
s.polynomialOrder = pol_deg;        %% Orden del polinomio
s.testRatio       = testratio;      %% Proporción de datos de prueba
data = Data(s);


s.networkParams.hiddenLayers    = hiddenLayers;
s.data                          = data;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;


opt = OptimizationProblem(s);
opt.solve();
%opt.plotRegressionResults(); 


if data.nFeatures == 2  % If you want to be asked for Features change it in "Data" Class
    opt.plotRegressionBoundary();
end

