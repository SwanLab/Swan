% This script generates multiple T data for various radii, applies SVD to
% extract the spatial modes U and generates V, which is then approximated
% via a NN in method B and via FE interpolation for method C.
close all;

radii = 0.5*(0:0.01:0.95);  % Radios para el conjunto de entrenamiento

% Crear generador de datos
dataGen = TrainingDataGenerator(radii);

fprintf('\nGenerando datos de entrenamiento\n');
dataGen.generateData(true);  % false = no calcular SVD (uses same method as for method a with direct NN)

dataGen.exportSVDToMAT('SVD_Results.mat','')
dataGen.exportSVDToCSV('');
% to regenerate T, apply T = U.*S*V'

%% Method B- SVD + NN to predict radial dependence

% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 40;
lambda          = 0.0;
learningRate    = 0.2;
hiddenLayers    = 128 .* ones(1, 6);

% INITIALIZATION - load model parameters
s.polynomialOrder = pol_deg;
s.testRatio       = testratio;
s.networkParams.hiddenLayers    = hiddenLayers;
s.optimizerParams.learningRate  = learningRate;
s.optimizerParams.maxEpochs     = 1000; % 1000 is the best option, but we use 100 to pass the tutorial quickly
s.costParams.lambda             = lambda;
s.costParams.costType           = 'L2';

s.networkParams.HUtype = 'ReLU';
s.networkParams.OUtype = 'linear';

s.fileName = 'VTrainingData.csv';
s.xFeatures = 1;
s.yFeatures = 2:(dataGen.getNumberOfModes()+1);

% Load Data
data   = Data(s);
s.data = data;

% Train the model
vPredictorNn = OptimizationProblemNN(s);
vPredictorNn.solve();
vPredictorNn.plotCostFnc();

save('VpredictorNN.mat', 'vPredictorNn'); % Save the model
Xin = [0.25];

Y = vPredictorNn.computeOutputValues(Xin);


%% Method C: SVD + HighOrder FE to predict radial dependence
degree = 20;  % High Order
VTrainingData = load('VTrainingData.csv');
r = VTrainingData(:,1);
v1 = VTrainingData(:,2);
v2 = VTrainingData(:,3);
p1 = polyfit(r, v1, degree);
p2 = polyfit(r, v2, degree); %not the best approach, quite some error as badly conditioned

% Save only these two tiny vectors (21 coefficients each)
save('PolynomialCoeffs.mat', 'p1', 'p2', 'degree');

% Online:
r_new = 0.3; 
v_new = [polyval(p1, r_new); polyval(p2, r_new)];

