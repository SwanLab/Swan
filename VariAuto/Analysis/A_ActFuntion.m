clear;clc;close all;addpath ../Codes;

%% INITIALIZATION
% Data
data          = Data('../Datasets/MNIST.csv',30,1);
% Network
hiddenlayers  = [500,150];
structure = [data.nFeatures,hiddenlayers,data.nLabels];
lambda        = 0;
% Trainer
learningRate      = 0.01;
opt.optTolerance  = 1*10^-10;
opt.maxevals      = 5000;
opt.maxepochs     = 5000;
opt.earlyStop     = 10;
opt.time          = Inf([1,1]);
opt.fv            = 10^-4;

%% HIDDEN UNIT ANALYSIS
actFCN = ["sigmoid","tanh","ReLU"];
l = 1;
historyFV   = zeros([l,length(actFCN)]);
historyTE   = zeros([l,length(actFCN)]);
historyTIME = zeros([l,length(actFCN)]);

for i = 1:length(actFCN)
    network   = Network(data,structure,'-loglikelihood',actFCN(i),'softmax',lambda);
    optimizer = Trainer.create(network,'SGD',learningRate,0,200,opt,'static');
    for j = 1:l
        network.computeInitialTheta();
        optimizer.train();

        historyTIME(j,i) = toc;
        [~,y_pred]       = max(network.getOutput(data.Xtest),[],2);
        [~,y_target]     = max(data.Ytest,[],2);
        testError        = mean(y_pred ~= y_target);
        historyFV(j,i)   = network.cost;
        historyTE(j,i)   = testError;
        fprintf('%d %f %f %f \n',j,historyFV(j,i),historyTE(j,i),historyTIME(j,i))
    end
    FV      = round(mean(historyFV,1),4);
    FV_sd   = std(historyFV,0,1);
    TE      = round(mean(historyTE,1),4);
    TE_sd   = std(historyTE,0,1);
    TIME    = round(mean(historyTIME,1),4);
    TIME_sd = std(historyTIME,0,1);
end

%% OUTPUT UNIT ANALYSIS
outFCN = ["sigmoid","softmax"];
l = 1;
historyFV2   = zeros([l,length(outFCN)]);
historyTE2   = zeros([l,length(outFCN)]);
historyTIME2 = zeros([l,length(outFCN)]);

for i = 1:length(outFCN)
    network   = Network(data,structure,'-loglikelihood-softmax','ReLU',outFCN(i),lambda);
    optimizer = Trainer.create(network,'SGD',learningRate,0,200,opt,'static');
    for j = 1:l
        network.computeInitialTheta();
        optimizer.train();

        historyTIME2(j,i) = toc;
        [~,y_pred]        = max(network.getOutput(data.Xtest),[],2);
        [~,y_target]      = max(data.Ytest,[],2);
        testError         = mean(y_pred ~= y_target);
        historyFV2(j,i)   = network.cost;
        historyTE2(j,i)   = testError;
        fprintf('%d %f %f %f \n',j,historyFV2(j,i),historyTE2(j,i),historyTIME2(j,i))
    end
    FV2      = round(mean(historyFV2,1),4);
    FV2_sd   = std(historyFV2,0,1);
    TE2      = round(mean(historyTE2,1),4);
    TE2_sd   = std(historyTE2,0,1);
    TIME2    = round(mean(historyTIME2,1),4);
    TIME2_sd = std(historyTIME2,0,1);
end