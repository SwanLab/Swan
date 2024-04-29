clear;clc;close all;addpath ../Codes;

%% INITIALIZATION
% Data choose between 32x32 or dct
s.fileName = '../Datasets/MNIST.csv';
s.polynomialOrder = 1;
s.testRatio       = 30;
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
