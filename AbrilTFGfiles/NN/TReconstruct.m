% Reconstruct T with a NN

clc 
clear


%% Loading Variables
% Load T NN
filename1='T_NN_test.mat';
filePath1 = fullfile('AbrilTFGfiles', 'NN', filename1);
load(filePath1);

% Mesh load
filename2='UL_r0_0000-20x20.mat';
filePath2 = fullfile('AbrilTFGfiles', 'DataVariables', filename2);
load(filePath2, 'mesh');

%% T reconstruct

r=0.5;
T=zeros(mesh.nnodes*mesh.ndim,8);

for j=1:8                       % Constructs the 8 columns    
    Taux2=[];
    for i=1:size(mesh.coord,1)  % Evaluates all the coordenates
        dataInput=[r,mesh.coord(i,:)];  %
        Taux1=T_NN(1,j).computeOutputValues(dataInput).';
        Taux2=cat(1,Taux2,Taux1);
    end
    T(:,j)=Taux2;
end

T2=T;



