% Reconstruct T with a NN and compare if it is correctly assemblied

clc 
clear


%% Loading Variables
% Load T NN
filename1='T_NN.mat';
filePath1 = fullfile('AbrilTFGfiles', 'NN', filename1);
load(filePath1,'T_NN');

% load the T real and mesh
fileName3="UL_r0_0000-50x50.mat";
filePath1 = fullfile('AbrilTFGfiles', 'Data','50x50', fileName3);
load(filePath1,"T","R","mesh");



%% T reconstruct
T_trained=zeros(mesh.nnodes*mesh.ndim,8);


for j=1:8 % Constructs the 8 columns    
    Taux2=[];
    for i=1:size(mesh.coord,1)  % Evaluates all the coordenates
        dataInput=[R,mesh.coord(i,:)];  
        Taux1=T_NN{1,j}.computeOutputValues(dataInput).';
        Taux2=cat(1,Taux2,Taux1);
    end
    T_trained(:,j)=Taux2;
end


%% COMPARISON ANALYTICAL VS TRAINED

diff=abs(T_trained-T);
ErrMax=max(max(diff));
disp('Error');
disp(ErrMax);


%% Results Exportation to paraview file
z.mesh      = mesh;
z.order     = 'P1';


for i=1:8
  z.fValues   = reshape(Tprova(:,i),[mesh.ndim,mesh.nnodes])';
  uFeFun = LagrangianFunction(z);%
  fileName = ['r03_T_NN' num2str(i)];
  centroids=computeCentroid(mesh);
  CoarsePlotSolution(uFeFun, mesh,[],fileName, R, centroids);
  %uFeFun.print(fileName,'Paraview');
end


for i=1:8
  z.fValues   = reshape(T_trained(:,i),[mesh.ndim,mesh.nnodes])';
  uFeFun = LagrangianFunction(z);%
  fileName = ['r03_T_NN' num2str(i)];
  centroids=computeCentroid(mesh);
  CoarsePlotSolution(uFeFun, mesh,[],fileName, R, centroids);
  %uFeFun.print(fileName,'Paraview');
end


function [centroids] =computeCentroid(mesh)
    x0=mean(mesh.coord(:,1));
    y0=mean(mesh.coord(:,2));
    centroids = [x0,y0];
end