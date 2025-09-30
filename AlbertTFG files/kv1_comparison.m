clc;
clear;
close all;


filePathData = fullfile('AlbertTFG files/mat files/Full/csvFiles_case1/K_r0.1000.csv');
NNname   = fullfile(['AlbertTFG files/opt Kv1 test5 learning0_01 hidden256x1_2 error1e-5.mat']);

data = readmatrix(filePathData);

opt = load(NNname);
opt = opt.opt;


error = [];

for i = 1:size(data,1)
    predicted = opt.computeOutputValues(data(i,1:8));
    error = cat(1, error, data(i,9:end)-predicted );
end

bar3c(abs(error));
zlabel("abs(diff)")
title("Abs of difference between real and predicted")
colormap('turbo');
clim([min(min(abs(error))), max(max(abs(error)))])
colorbar()

xlabel("Position i in array")
ylabel("Position j in array")