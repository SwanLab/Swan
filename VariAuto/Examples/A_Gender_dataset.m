%% Testing with a sample nn
clc;
clear;
close all;

%% Initialization of hyperparameters
pol_deg         = 1;
testratio       = 30;  
lambda          = 0;
learningRate    = 0.1;
momentum        = 0.9;
batch           = 200;
hiddenlayers    = [10,15];

%% Loading of files/datasets
datasets = load("datasets.mat").datasets1;
disp('Datsets available:')
for i = 1:length(datasets)
    fprintf('%d - %s \n',i,datasets(i))
end
fileN = 'Gender.csv';%datasets(input('Choose: '));

s.features = 1:2;
s.fileName        = fileN;
s.testRatio       = testratio;
s.polynomialOrder = pol_deg;
data  = Data(s);

%% Overview of the data (divided in men and women)
% overv = load(fileN);
% hold on
% for i = 1:size(overv,1)
%     if overv(i,3) == 1
%         plot(overv(i,1),overv(i,2),'og') % Men
%     elseif overv(i,3) == 2 
%         plot(overv(i,1),overv(i,2),'or') % Women
%     end
% end
% hold off
% xlabel("Height (cm)")
% ylabel("Weight (kg)")

%% Create Network and trainer Objects
structure = [data.nFeatures,hiddenlayers,data.nLabels];
% network   = Network(data,structure);
% network = Network(data,structure,'-loglikelihood','ReLU','softmax',lambda);

%% Run Optimization Problem
p.data            = data;
p.structure       = structure;
p.optimizerParams.learningRate = learningRate;
p.costParams.lambda = lambda;
p.networkParams.hiddenLayers = hiddenlayers;

optProblem   = OptimizationProblem(p);
% opt.optTolerance  = 1*10^-8; opt.maxevals      = 100;
% optProblem.maxepochs     = 100; 
% opt.earlyStop     = 10;
% opt.time          = Inf([1,1]); opt.fv         = 10^-4;
% nplt              = 1;
% optimizer       = Trainer.create(network,'SGD',learningRate,momentum,batch,opt,'static',nplt);

optProblem.solve();
%% RUN & Possible functions
optProblem.plotCostFnc();
% optProblem.plotConections(); NOT WORKING
% optProblem.plotBoundary('contour'); NOT WORKING
% optProblem.plotSurface(); NO ENTIENDO UTILIDAD
optProblem.plotConfusionMatrix();