% Reconstruct K with a NN and compare if it is correctly assemblied

clc;
clear;
close all;

%% Case parameters
p.Sampling   ='Isolated';     %'Isolated'/'Oversampling'
p.Inclusion  ='Material';    %'Material'/'Hole'/'HoleRaul
p.nelem      = 50;
meshName    =  p.nelem+"x"+p.nelem;


%% LOAD DATA

% load the trained NN 
filePath = fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,"K_NN.mat");
load(filePath1);

% load K coarse dataset
fileName2="DataK.csv";
filePath2 = fullfile('AbrilTFGfiles',fileName2);
load(filePath2);

% load the K real
fileName3="UL_r0_5000-20x20.mat";
filePath1 = fullfile('AbrilTFGfiles', 'Data', '20x20',fileName3);
load(filePath1,"K","R");



%% RESTRUCTURE DATA

K_aux=K_NN.computeOutputValues(R);
n=8;
K_NN=zeros(n);

idx=1
for i=1:n
    for j=i:n
        K_NN(i,j)=K_aux(idx);
        idx=idx+1;
    end
end

K_NN=K_NN+triu(K_NN,1).';

%% COMPARISON ANALYTICAL VS TRAINED

diff=abs(K_NN-K)
ErrMax=max(max(diff));
disp('Error');
disp(ErrMax);



