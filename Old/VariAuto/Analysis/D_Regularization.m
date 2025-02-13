clear;clc;close all;addpath ../Codes;

%% INITIALIZATION
% Data
data          = Data('../Datasets/Microchip.csv',30,10);
% Network
hiddenlayers = [];
structure    = [data.nFeatures,hiddenlayers,data.nLabels];
% Trainer
learningRate      = 0.1;
opt.optTolerance  = 1*10^-4;
opt.maxevals      = 5000;
opt.maxepochs     = 30;
opt.time          = Inf([1,1]);
opt.fv            = 0;

tr     = [25,50,75];
n = length(tr);
l = 10;
%% L2 ANALYSIS
lambda = [10^-8,10^-6,10^-4,10^-2,10^0];
m = length(lambda);
fv = zeros([l,n]);
te = zeros([l,n]);

for k = 1:m
    network   = Network(data,structure,'-loglikelihood','ReLU','softmax',lambda(k));
    optimizer = Trainer.create(network,'SGD',learningRate,0);
    for i = 1:n
        for j = 1:l
            network.computeInitialTheta();
            optimizer.train();
            fprintf('%d',i)

            [~,y_pred] = max(network.getOutput(data.Xtest),[],2);
            [~,y_target] = max(data.Ytest,[],2);
            fv(j,i,k) = network.cost;
            te(j,i,k) = mean(y_pred ~= y_target);
        end
        TE(k,i) = mean(te(:,i,k));
        TE_sd(k,i) = std(te(:,i,k));
    end
end

%% EARLY STOP ANALYSIS
ES     = [1 5 10 50 100 150 200 250];
m = length(ES);
fv = zeros([l,n]);
te = zeros([l,n]);
for k = 1:m
    network   = Network(data,structure);
    opt.earlyStop = ES(k);
    optimizer = Trainer.create(network,'SGD',learningRate,0,200,opt,'static');
    for i = 1:n
        for j = 1:l
            network.computeInitialTheta();
            optimizer.train();
            fprintf('%d',i)

            [~,y_pred] = max(network.getOutput(data.Xtest),[],2);
            [~,y_target] = max(data.Ytest,[],2);
            fv(j,i,k) = network.cost;
            te(j,i,k) = mean(y_pred ~= y_target);
        end
        TE2(k,i) = mean(te(:,i,k));
        TE2_sd(k,i) = std(te(:,i,k));
    end
end

for k = 1:length(tr)
    figure(1)
    errorbar(lambda,TE(:,k),TE_sd(:,k))
    legend
    hold on

    figure(2)
    errorbar(ES,TE2(:,k),TE2_sd(:,k))
    legend
    hold on
end