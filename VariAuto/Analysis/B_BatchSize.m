clear;clc;close all;addpath ../Codes;

%% INITIALIZATION
% Data
data          = Data('../Datasets/MNIST.csv',30,1);
% Network
hiddenlayers = [500,150];
structure    = [data.nFeatures,hiddenlayers,data.nLabels];
lambda       = 0;
network      = Network(data,structure);
% Trainer
learningRate      = 0.01;
momentum          = [0,0.9,0.9];
opt.optTolerance  = 1*10^-10;
opt.maxevals      = 5000;
opt.maxepochs     = 5000;
opt.earlyStop     = 5000;
% These two for target time
opt.time          = 60;
opt.fv            = 0;
% These two for target function value
% opt.time          = Inf([1,1]);  
% opt.fv            = 0.01;

%% BATCH ANALYSIS
optType     = ["SGD","Nesterov","RMSProp"];
batch       = [1 5 10 50 100 200 500 1000 2000 4000 8000];
l           = 10;
historyFV   = zeros([l,length(batch),length(optType)]);
historyTE   = zeros([l,length(batch),length(optType)]);
historyTIME = zeros([l,length(batch),length(optType)]);

for k = 1:length(optType)
    for i = 1:length(batch)
        optimizer = Trainer.create(network,optType(k),learningRate,momentum(k),batch(i),opt,'static');
        for j = 1:l
        network.computeInitialTheta();
        optimizer.train();
        
        historyTIME(j,i,k) = toc;
        [~,y_pred]         = max(network.getOutput(data.Xtest),[],2);
        [~,y_target]       = max(data.Ytest,[],2);
        testError          = mean(y_pred ~= y_target);
        historyFV(j,i,k)   = network.cost;
        historyTE(j,i,k)   = testError;
        fprintf('%d %d TE %f\n',i,j,historyTE(j,i))
        if j >= 3
            if mean(historyTIME(1:j,i,k), 1) > 950
                historyTIME(j+1:end,i,k) = 1000;
                break
            end
        end
       end
    end
end
FV      = round(mean(historyFV,1),4);
FV_sd   = std(historyFV,0,1);
TE      = round(mean(historyTE,1),4);
TE_sd   = std(historyTE,0,1);
TIME    = round(mean(historyTIME,1),4);
TIME_sd = std(historyTIME,0,1);

for k = 1:length(optType)
    figure(1)
    errorbar(batch,FV(:,:,k),FV_sd(:,:,k))
    figure(2)
    errorbar(batch,TE(:,:,k),TE_sd(:,:,k))
    figure(3)
    errorbar(batch,TIME(:,:,k),TIME_sd(:,:,k))
    legend("SGD","Nesterov","RMSProp")    
end



