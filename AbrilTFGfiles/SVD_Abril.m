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
nBasis=10;

basis = U(:,1:nBasis);
sValues = S(1:nBasis,1:nBasis);
q = V(:,1:nBasis)*sValues;

table=[r q];
%QFileName = fullfile('AbrilTFGfiles', 'DataQ.csv');
%writematrix(table,QFileName);



%% Test NN
NNname=fullfile("AbrilTFGfiles/NN/Q_NN.mat");
load(NNname);

rad=0.3;
qNN=Q_NN.computeOutputValues(rad).';
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

tiledlayout(10,5,'TileSpacing','compact','Padding','compact');
for i=1:50
    ax=nexttile;
    plot(r,V(:,i), 'LineWidth', 1.5);
    xlabel('r');
    ylabel("V"+i);
    title("V"+ i+" ; r="+r(i,1));
    grid on
end

% Plot the S
figure
plot(log(diag(S)),'LineWidth',1.5);
title("S singular values");
ylabel('value');
xlabel('column');

% Plot U
%s.mesh=mesh;
%s.order='P1';
%n=10;
%
%Ufun=cell(1,n);
%for i=1:n
%    s.fValues=U(:,i);
%    Ufun{1,i}=LagrangianFunction(s);
%    Ufun{1,i}.plot();
%end

function centroids=computeCentroid(mesh)
    x0=mean(mesh.coord(:,1));
    y0=mean(mesh.coord(:,2));
    centroids = [x0,y0];
end
    

