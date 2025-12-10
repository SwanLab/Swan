%SVD FOR THE NN OPTION 3

clc;
clear;
close all;

%% Load Data
directory = './AbrilTFGfiles/Data/50x50';
files     = dir(fullfile(directory, 'UL_*.mat'));

r=zeros(length(files),1);
for k = 1:1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);
    data=load(filePath);
    if k==1
        T=zeros(data.mesh.nnodes*data.mesh.ndim*8,length(files));
    end
    T(:,k) = data.T(:);
    r(k,1)=data.R;
end

%% Data to train
[U,S,V]    = svd(T,"econ");
nBasis=20;

basis = U(:,1:nBasis);
sValues = S(1:nBasis,1:nBasis);
q = V(:,1:nBasis)*sValues;

table=[r q];

%fileName=fullfile("AbrilTFGfiles","SVD","SVD.mat");
%save(fileName,"U","S","V","r");
%QFileName = fullfile('AbrilTFGfiles', 'DataQ.csv');
%writematrix(table,QFileName);
%% Graphics

% PLot the V columns grouped in 10
%step=10;
%Nwindow=ceil(size(V,2)/step);
%idx=1;

%for j=1:Nwindow
%    figure('Position',[75 100 1400 600]);
%    tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
%    for i=1:step
%        ax=nexttile;
%        plot(r,V(:,idx), 'LineWidth', 1.5);
%        xlabel('r');
%        ylabel("V(:,"+idx);
%        title("V-"+ idx+" ; r="+r(1,idx));
%        grid on
%        idx=idx+1;
%    end
%end

tiledlayout(3,5,'TileSpacing','compact','Padding','compact');
for i=1:15
    ax=nexttile;
    plot(r,V(:,i), 'LineWidth', 1.5);
    xlabel('r');
    ylabel("V"+i);
    title("V"+ i);
    legend('Exact')
    grid on
end

% Plot the S
figure
plot(log(diag(S)),'LineWidth',1.5);
title("S singular values");
ylabel('value');
xlabel('column');


%% Export NN paraview file
NNname=fullfile("AbrilTFGfiles/NN/Q_NN.mat");
load(NNname);



rad=0.3;
qNN=zeros(15,1);
for i=1:size(Q_NN,2)
    qNN(i)=Q_NN{i}.computeOutputValues(rad).';
end

T_NN=basis*qNN;

u=reshape(T_NN,[],8);

% Export vtu file
 z.mesh      = data.mesh;
 z.order     = 'P1';
 for i=1:8
   z.fValues   = reshape(u(:,i),[data.mesh.ndim,data.mesh.nnodes])';
   uFeFun = LagrangianFunction(z);%
   fileName = strrep("r" + num2str(rad), '.', '_')+ "_SVD_NN_Training" +num2str(i);
   centroids=computeCentroid(data.mesh);
   CoarsePlotSolution(uFeFun, data.mesh,[],fileName, rad, centroids);
 end

function centroids=computeCentroid(mesh)
    x0=mean(mesh.coord(:,1));
    y0=mean(mesh.coord(:,2));
    centroids = [x0,y0];
end
    