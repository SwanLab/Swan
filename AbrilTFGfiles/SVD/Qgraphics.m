%% Q graphics Comparison

% 1. SVD Solution
% 2. HO FE aproximation
% 3. Hybrid SVD + NN

clc; clear;

%% LOAD DATA
p.nelem=20;
p.Training  = 'EIFEM';         % 'EIFEM'/'Multiscale'
p.Sampling ='Oversampling';     %'Isolated'/'Oversampling'
p.Inclusion='Material';    %'Material'/'Hole'/'HoleRaul
meshName    = p.nelem+"x"+p.nelem;

% NN
NNname=fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,meshName,"Q_NN.mat");
load(NNname);

% High Order function
HOname=fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,meshName,"HOfunction.mat");
load(HOname,"fR");

% SVD
SVDname=fullfile("AbrilTFGfiles","Data",p.Training,p.Inclusion,p.Sampling,meshName,"SVD.mat");
SVDdata=load(SVDname);


%% SHAPES DATA

r=SVDdata.r;

svdValues=SVDdata.V*SVDdata.S;

% NN
rFull = Data.buildModel(r,pol_deg);
NNvalues=Q_NN.computeOutputValues(rFull);


% HO (Esta a escalat al domini -1 1)
HOvalues=zeros(size(r,1),size(fR,2));
for i=1:size(r,1)
    for j=1:size(fR,2)
        HOvalues(i,j)=fR(j).evaluate(2*(r(i) - min(r)) / (max(r) - min(r)) - 1);
    end
end


%% Graphics

% Q graphic
pos=[271,316,1451,568];
tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
for i=1:10
    ax=nexttile;
    plot(r,svdValues(:,i),r,HOvalues(:,i), r,NNvalues(:,i),'LineWidth', 1.1);
    xlabel('p');
    ylabel("\alpha"+i);
    title("\alpha"+ i);
    legend('Dataset','HO','NN');
    grid on
end
set(gcf, 'Position', pos);