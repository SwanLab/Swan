clear;clc;close all;addpath ../Codes;

%% INITIALIZATION
% Data choose between 32x32 or dct
data          = Data('../Datasets/FlowersDCT.csv',30,1);
% data          = Data('../Datasets/Flowers3232.csv',30,1);
% Network
hiddenlayers  = [1000,500,150];
structure = [data.nFeatures,hiddenlayers,data.nLabels];
network   = Network(data,structure);
% Trainer
learningRate      = 0.01;
momentum          = 0.5;
batch             = 200;
opt.optTolerance  = 1*10^-10;
opt.maxevals      = 500;
opt.maxepochs     = 500;
opt.earlyStop     = 500;
opt.time          = Inf([1,1]);
opt.fv            = 10^-4;
optimizer         = Trainer.create(network,'Nesterov',learningRate,momentum,batch,opt,'static');

optimizer.train();
network.plotConfusionMatrix();
%% UNCOMMENT FOR PLOTTING SOME WRONG IMAGES
% [~,OUT]   = max(network.getOutput(data.Xtest),[],2);
% [~,TAR] = max(data.Ytest,[],2);
% err = TAR ~= OUT;
% 1 - sum(err)/length(err)
% IDX = [];
% for i = 1:size(data.Xtest,1)
%     if err(i) == 1
%         IDX = [IDX,i];
%     end
% end
% 
% amin = min([-13683.9464995775,-14685.1858341210,-14529.1420506253,-14297.2434420311,-13075.3659637298]);
% amax = max([42736.9111111111,45371.0833333333,44924.5555555555,42145.2055555555,45160.7611111111]);
% mins = min(data.data(:,1:end-1),[],1);
% maxs = max(data.data(:,1:end-1),[],1);
% X = data.Xtest.*(maxs-mins) + mins;
% for k = 1:10
%     reconstruct = zeros([180,180,3]);
%     reconstructIDCT = zeros([180,180,3]);
%     cont = 1;
%     for c = 1:3
%         for i = 1:45
%             for j = 1:45
%                 if i + j <=46
%                     reconstruct(j,i,c) = X(IDX(k),cont)*(amax-amin) + amin;
%                     cont = cont+1;
%                 end
%             end
%         end
%         reconstructIDCT(:,:,c) = idct2(reconstruct(:,:,c));
%     end
%     figure(k)
%     imshow(rescale(reconstructIDCT))
%     tlt = strcat('Predicted = ',num2str(OUT(IDX(k))),' , Target = ',num2str(TAR(IDX(k))));
%     title(tlt)
% end