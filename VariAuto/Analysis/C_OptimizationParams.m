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
opt.optTolerance  = 1*10^-10;
opt.maxevals      = 5000;
opt.maxepochs     = 30;
opt.earlyStop     = 5000;
opt.time          = Inf([1,1]);
opt.fv            = 0;
nPlot             = 1;   

%% OPTPARAMS ANALYSIS
optType           = ["SGD","Nesterov","RMSProp"];
learningRate = [0.0001 0.001 0.01 0.1 1; 
                0.01 0.01 0.01 0.01 0.01;
                0.01 0.01 0.01 0.01 0.01];
momentum     = [0  0   0   0; 
                0. 0.5 0.9 0.99;
                0. 0.5 0.9 0.99];
n = [size(learningRate,2),size(momentum,2),size(momentum,2)];
y = cell([1,length(optType)]);
x = cell([1,length(optType)]);
for i = 1:length(optType)
    y{i} = zeros([n(i),opt.maxepochs-1]);
    x{i} = zeros([n(i),opt.maxepochs-1]);
end
for k = 1:length(optType)
    for i = 1:n(k)
        network.computeInitialTheta();
        optimizer = Trainer.create(network,optType(k),learningRate(k,i),momentum(k,i),200,opt,'static',nPlot);
        optimizer.train();
        fprintf('%d',i)

        dataObjsY = findobj(figure((i-1)*2+1),'-property','YData');
        y{k}(i,:) = dataObjsY(1).YData;
        dataObjsX = findobj(figure((i-1)*2+1),'-property','XData');
        x{k}(i,:) = dataObjsX(1).XData;
        figure(k+10)
        plot(x{k}(i,:),y{k}(i,:),'d-')
        hold on
    end
    hold off
end
