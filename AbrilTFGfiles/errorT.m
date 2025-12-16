%% T error Graphics Comparison

% This script has the purpose to visualize the error between the different
% options along the radius. It is initially intended just to study for the
% mesh 20x20 case

% 1. T NN 
% 2. HO FE aproximation
% 3. Hybrid SVD + NN

clc; clear;

%% LOAD DATA

% 1. NN
filename='T_NN3.mat';
filePath = fullfile('AbrilTFGfiles', 'NN', filename);
load(filePath,'T_NN');

% 2. High Order function
HOname=fullfile("AbrilTFGfiles/SVD/HOfunction.mat");
load(HOname,"fT","deim");

% 3. SVD +NN
NNname=fullfile("AbrilTFGfiles/NN/Q_NN2.mat");
load(NNname);
U=deim.basis(:,1:10);


% Dataset
directory = './AbrilTFGfiles/Data/50x50';
files = dir(fullfile(directory, 'UL_*.mat'));
i=1;
for k = 1:1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    data=load(filePath);
    T(:,i) = data.T(:);  % This stores each file's contents in the cell array 'allData'
    %disp(['Loaded: ', files(k).name]);  % Display the file being loaded
    i=i+1;
end

mesh=data.mesh;

%% RECONSTRUCT T DATA
r=0:0.02:0.98;

T1= zeros(mesh.nnodes*mesh.ndim*8,length(r));
T2= zeros(mesh.nnodes*mesh.ndim*8,length(r));
T3= zeros(mesh.nnodes*mesh.ndim*8,length(r));

Tprova=computeT2(mesh,0.3,T_NN);
% 1. NN
for i=1:length(r)
    aux=computeT2(mesh,r(i),T_NN);
    T1(:,i)=aux(:);
end

% 2. High Order function
for i=1:length(r)
    aux=fT(r(i));
    T2(:,i)=aux(:);
end

% 3. SVD + NN
for i=1:length(r)
    rFull = Data.buildModel(r(i),6);
    q=Q_NN.computeOutputValues(rFull).';
    T3(:,i)=U*q;
end

%% CALCULATE ERROR

err1=vecnorm(abs(T-T1))/norm(T);
err2=vecnorm(abs(T-T2))/norm(T);
err3=vecnorm(abs(T-T3))/norm(T);


%% PLOT ERROR

figure
plot(r,err1,r,err2,r,err3,LineWidth=1.5);
legend("NN","HO","SVD+NN");
title("Error vs r for 50x50 mesh");
xlabel('r');
ylabel('error');






%% Functions
function T_trained=computeT(mesh,R,T_NN)
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
end

function T_trained=computeT2(mesh,R,T_NN)
    T_trained=[];
        Taux2=[];
        for i=1:size(mesh.coord,1)  % Evaluates all the coordenates
            dataInput=[R,mesh.coord(i,:)];
            dataFull=Data.buildModel(dataInput,6);
            Taux1=T_NN.computeOutputValues(dataFull).';
            Taux2=reshape(Taux1,2,[]);
            T_trained=[T_trained;Taux2];
        end
end
