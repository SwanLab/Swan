% Reconstruct T with a NN and compare if it is correctly assemblied

clc 
clear


%% Loading Variables
% Load T NN
filename1='T_NN.mat';
filePath1 = fullfile('AbrilTFGfiles', 'NN', filename1);
load(filePath1,'T_NN');

% load the T real and mesh
fileName3="UL_r0_5000-50x50.mat";
filePath1 = fullfile('AbrilTFGfiles', 'Data','50x50', fileName3);
load(filePath1,"T","R","mesh");


%% T reconstruct

T_trained=[];
    Taux2=[];
    for i=1:size(mesh.coord,1)  % Evaluates all the coordenates
        dataInput=[R,mesh.coord(i,:)];
        dataFull=Data.buildModel(dataInput,6);
        Taux1=T_NN.computeOutputValues(dataFull).';
        Taux2=reshape(Taux1,2,[]);
        T_trained=[T_trained;Taux2];
    end

exportT_weakInclusion(T_trained,R,mesh,"R05testFunction");

%% COMPARISON ANALYTICAL VS TRAINED

diff=abs(T_trained-T);
ErrMax=max(max(diff));
disp('Error');
disp(ErrMax);
