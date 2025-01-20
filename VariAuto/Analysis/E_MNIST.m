clear;clc;close all

%% INITIALIZATION
% Data choose between 32x32 or dct
s.fileName = 'MNIST.csv';
s.polynomialOrder = 1;
s.testRatio       = 30;
s.features        = 1:784;
data = Data(s);


hiddenLayers  = [500,100,500];

learningRate      = 0.04;
lambda = 0;

s.networkParams.hiddenLayers    = hiddenLayers;
s.data                          = data;
s.optimizerParams.learningRate  = learningRate;
s.costParams.lambda             = lambda;

opt = OptimizationProblem(s);
opt.solve();

opt.plotImage(2001);
opt.plotCostFnc
