%% PRECONDITIONING SCRIPT

%This script is intented to execute the coarseTesting_Abril for varius
%cases and create plots that summarise all the results in a better way

%% MAIN INPUTS

% case parameters
s.Training  = 'EIFEM';            % 'EIFEM'/'Multiscale'
s.Inclusion = 'Material';         % 'Hole'/'Material'/'HoleRaul'   --> Hole: just for constant r
s.Sampling  = 'Oversampling';     % 'Isolated'/'Oversampling'
s.Option    = [];                 % 'Dataset'/'NN'/'HO'/ 'Hybrid'
s.nelem     =  20;                %  Mesh refining
s.Print     = false;

% UNIFORM DISTRIBUTION
s.r = ones(3,10)*0.3;
% s.r = [0.1,0.2,0.3,0.4,0.5
%          0.1,0.2,0.3,0.4,0.5
%          0.1,0.2,0.3,0.4,0.5];

%NON-UNIFORM DISTRIBUTION
% s.r=[0.1,0.2,0.3
%      0.4,0.5,0.6
%      0.7,0.8,0.8];


%% DATASET OPTION
s.Option = 'Dataset';
DSet= CoarseTesting_2D(s);
DSet.compute();

%% HO OPTION
s.Option = 'HO';
HO= CoarseTesting_2D(s);
HO.compute();

%% HYBRID OPTION
s.Option = 'Hybrid';
Hybd= CoarseTesting_2D(s);
Hybd.compute();

%% NN OPTION
s.Option = 'NN';
NN= CoarseTesting_2D(s);
NN.compute();

%% COMPARE PLOTS

% RESIDUAL
figure
plot(DSet.residualPCG,'linewidth',1.5)
hold on
plot(HO.residualPCG,'linewidth',1.5)
hold on
plot(Hybd.residualPCG,'linewidth',1.5)
hold ON
plot(NN.residualPCG,'linewidth',1.5)
set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual Evolution EIFEM Oversampling")
legend({'CG Dataset', 'CG HO', 'CG Hybrid','CG NN'});