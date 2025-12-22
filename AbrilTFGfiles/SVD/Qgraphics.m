%% Q graphics Comparison

% 1. SVD Solution
% 2. HO FE aproximation
% 3. Hybrid SVD + NN

clc; clear;

%% LOAD DATA
p.nelem=20;
p.Sampling ='Isolated';     %'Isolated'/'Oversampling'
p.Inclusion='Material';    %'Material'/'Hole'/'HoleRaul
p.nelem=50;
meshName    = p.nelem+"x"+p.nelem;

% NN
NNname=fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,"Q_NN.mat");
load(NNname);

% High Order function
HOname=fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,"HOfunction.mat");
load(HOname,"fR");

% SVD
SVDname=fullfile("AbrilTFGfiles","Data",p.Inclusion,p.Sampling,"SVD.mat");
SVDdata=load(SVDname);


%% SHAPES DATA

r=SVDdata.r;

svdValues=SVDdata.V*SVDdata.S;

% NN
rFull = Data.buildModel(r,6);
NNvalues=Q_NN.computeOutputValues(rFull);



%NNvalues=zeros(size(r,1),size(Q_NN,2));
%for j=1:size(Q_NN,2)     
%    rFul = Data.buildModel(r,9);
%    NNvalues(:,j)=Q_NN{j}.computeOutputValues(rFul);
%end


% HO (Esta a escalat al domini -1 1)
HOvalues=zeros(size(r,1),size(fR,2));
for i=1:size(r,1)
    for j=1:size(fR,2)
        HOvalues(i,j)=fR(j).evaluate(2*(r(i) - min(r)) / (max(r) - min(r)) - 1);
    end
end


%% Graphics

% Q graphic
tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
for i=1:10
    ax=nexttile;
    plot(r,svdValues(:,i),r,HOvalues(:,i), r,NNvalues(:,i),'LineWidth', 1);
    xlabel('r');
    ylabel("Q"+i);
    title("Q"+ i);
    legend('Exact','HO','NN');
    grid on
end