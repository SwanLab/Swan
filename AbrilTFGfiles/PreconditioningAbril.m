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

% Definition of Subdomains
% s.r= ones(10,10)*0.1;
% s.r = ones(5,15)*0.1;
% s.r = [0.1,0.2,0.3,0.4,0.5
%          0.1,0.2,0.3,0.4,0.5
%          0.1,0.2,0.3,0.4,0.5];
s.r=[0.1,0.2,0.3
     0.4,0.5,0.6
     0.7,0.8,0.9];


%% DATASET OPTION
s.Option = 'Dataset';
DSet= CoarseTesting_AbrilV2(s);
DSet.compute();
% DSet.PlotSolution();

%% HO OPTION
s.Option = 'HO';
HO= CoarseTesting_Abril(s);
HO.compute();
% HO.PlotSolution;

%% HYBRID OPTION
s.Option = 'Hybrid';
Hybd= CoarseTesting_Abril(s);
Hybd.compute();
% Hybd.PlotSolution();

%% NN OPTION
s.Option = 'NN';
NN= CoarseTesting_Abril(s);
NN.compute();
% NN.PlotSolution();

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
% hold on
% plot(DSet.residualCG,'linewidth',1.5)
set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual 10x10, r=0.1 - EIFEM Oversampling")
legend({'PCG Dataset', 'PCG HO', 'PCG Hybrid','PCG NN'});
% legend({'CG Dataset', 'CG HO', 'CG Hybrid','CG NN','CG + ILU-EIFEM-ILU'});