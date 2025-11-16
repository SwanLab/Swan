% Reconstruct T with a NN and compare if it is correctly assemblied

clc 
clear


%% Loading Variables
% Load T NN
filename1='T_NN_test.mat';
filePath1 = fullfile('AbrilTFGfiles', 'NN', filename1);
load(filePath1);

% load the T real and mesh
fileName3="UL_r0_5000-20x20.mat";
filePath1 = fullfile('AbrilTFGfiles', 'DataVariables', fileName3);
load(filePath1,"T","R","mesh");
%% T reconstruct


T_trained=zeros(mesh.nnodes*mesh.ndim,8);

for j=1:8 % Constructs the 8 columns    
    Taux2=[];
    for i=1:size(mesh.coord,1)  % Evaluates all the coordenates
        dataInput=[R,mesh.coord(i,:)];  
        Taux1=T_NN(1,j).computeOutputValues(dataInput).';
        Taux2=cat(1,Taux2,Taux1);
    end
    T_trained(:,j)=Taux2;
end

T2=T;

%% COMPARISON ANALYTICAL VS TRAINED

diff=abs(T_trained-T)
ErrMax=max(max(diff));
disp('Error');
disp(ErrMax);



