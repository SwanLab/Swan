%% Testing with a sample nn
clc;
clear;
close all;
%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 20;
lambda          = 0;
learningRate    = 0.05;
momentum        = 0.9; % Not used
batch           = 200; % Not used
hiddenlayers    = [];


%% Loading of files/datasets
fileN = 'finalEllipseDataset.csv';
s.features = 1:3;
s.fileName        = fileN;
s.testRatio       = testratio;
s.polynomialOrder = pol_deg;
data  = Data(s);

data.Ytest = data.Ytest(~isnan(data.Ytest(:,1)),:);
data.Xtest = data.Xtest(~isnan(data.Ytest(:,1)),:);
data.Xtrain = data.Xtrain(~isnan(data.Ytrain(:,1)),:);
data.Ytrain = data.Ytrain(~isnan(data.Ytrain(:,1)),:);
data.Ntest = size(data.Ytest,1);

%% Create Network and trainer Objects
structure = [size(data.Xtrain,2),hiddenlayers,data.nLabels];

%% Run Optimization Problem
p.data            = data;
p.structure       = structure;
p.optimizerParams.learningRate = learningRate;
p.costParams.lambda = lambda;
p.networkParams.hiddenLayers = hiddenlayers;
p.networkParams.costType     = 'L2';
p.networkParams.HUtype       = 'tanh';
p.networkParams.OUtype       = 'None';
optProblem   = OptimizationProblem(p);

tic;
optProblem.solve();
t = toc;

%% Additional plots
optProblem.plotCostFnc();

%% Error
% Evaluation for the train data
errTrain = optProblem.computeError(data.Xtrain, data.Ytrain);
% Evaluation for the test data
errTest = optProblem.computeError(data.Xtest, data.Ytest);

% %% Surfaces differences
X = data.Xtrain(:,1);
Y = data.Xtrain(:,2);
z_theo = data.Ytrain;
z_calc = optProblem.computeOutputValues(data.Xtrain);

% Define la cuadrícula en la que se interpolan los puntos
[Xq,Yq] = meshgrid(min(X):0.01:max(X), min(Y):0.01:max(Y));

z_vec_label = {'C_{11}','C_{12}','C_{13}','C_{22}','C_{23}','C_{33}'};
for i = 1:6
    Zq = griddata(X, Y, z_calc(:,i), Xq, Yq);
    Zp = griddata(X, Y, z_theo(:,i), Xq, Yq);
    
    % Crea el gráfico de superficie
    figure ();
    surf(Xq, Yq, Zq,'FaceColor','red');
    hold on
    surf(Xq, Yq, Zp,'FaceColor','green');
    xlabel('a')
    ylabel('b')
    zlabel(z_vec_label{i})
    lgd = legend('NN Calculation','Comp. Homog. Calculation');
    lgd.FontSize = 17;
end












