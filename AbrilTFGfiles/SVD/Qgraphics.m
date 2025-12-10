%% Q graphics Comparison

% 1. SVD Solution
% 2. HO FE aproximation
% 3. Hybrid SVD + NN

clc; clear;

%% LOAD DATA

% NN
NNname=fullfile("AbrilTFGfiles/NN/Q_NN.mat");
load(NNname);

% High Order function
HOname=fullfile("AbrilTFGfiles/SVD/HOfunction.mat");
load(HOname);

% SVD
SVDname=fullfile("AbrilTFGfiles/SVD/SVD.mat");
SVDdata=load(SVDname);


%% SHAPES DATA

r=SVDdata.r;

svdValues=SVDdata.V*SVDdata.S;

% NN
%NNvalues=zeros(size(r,1),20);
%for i=1:size(r,1)
%    NNvalues(i,:)=Q_NN.computeOutputValues(r(i))
%end


NNvalues=zeros(size(r,1),size(Q_NN,2));
for i=1:size(r,1)
    for j=1:size(Q_NN,2)
        NNvalues(i,j)=Q_NN{j}.computeOutputValues(r(i));
    end
end



% HO (Esta a escalat al domini -1 1)
HOvalues=zeros(size(r,1),size(fR,2));
for i=1:size(r,1)
    for j=1:size(fR,2)
        HOvalues(i,j)=fR(j).evaluate(2*(r(i) - min(r)) / (max(r) - min(r)) - 1);
    end
end


%% Graphics

% Q graphic
tiledlayout(3,5,'TileSpacing','compact','Padding','compact');
for i=1:15
    ax=nexttile;
    plot(r,svdValues(:,i),r,HOvalues(:,i), r,NNvalues(:,i),'LineWidth', 1);
    xlabel('r');
    ylabel("Q"+i);
    title("Q"+ i);
    legend('Exact','HO','NN');
    grid on
end