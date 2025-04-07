
clear;
clc;

file = 'E_AoA0_mpt.txt';
data = readmatrix(file);
realData = data(:,[1:4,end]);

sN.data.nFeatures    = 4;
sN.data.nLabels      = 1;
sN.hiddenLayers = 2;

n = Network(sN);
learningRate = 0.01;
momentum     = 0;

sT.type = 'SGD';
t = Trainer.create(sT);