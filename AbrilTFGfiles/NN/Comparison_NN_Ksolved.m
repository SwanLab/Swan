% COMPARISON OF THE RESULTS OF THE K_NN TRAINED VS NORMAL

clc;
clear;
close all;

%% LOAD DATA

%load the trained NN 
fileName1="K_NN.mat";
filePath1 = fullfile('AbrilTFGfiles', 'NN', fileName1);
load(filePath1);

%%load K coarse dataset
fileName2="Kdata.csv";
filePath2 = fullfile('AbrilTFGfiles',fileName2);
load(filePath2);


%% RESTRUCTURE DATA
K_auxiliar=opt.computeOutputValues(R);
n=8;
K_NN=zeros(n);

idx=1
for i=1:n
    for j=i:n
        K_NN(i,j)=K_auxiliar(idx);
        idx=idx+1;
    end
end

K_NN=K_NN+triu(K_NN,1).';

%% COMPARISON ANALYTICAL VS TRAINED
Err=max(max(abs(K_NN-K)));
disp('Error');
disp(Err);

%disp(abs(K_NN-K));


