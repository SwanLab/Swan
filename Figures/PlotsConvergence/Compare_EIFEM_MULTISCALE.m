%% PRECONDITIONING SCRIPT

%This script is intented to execute the coarseTesting_Abril for varius
%cases and create plots that summarise all the results in a better way

%% MAIN INPUTS

% case parameters
s.Training  = [];                 % 'EIFEM'/'Multiscale'
s.Inclusion = 'Material';         % 'Hole'/'Material'/'HoleRaul'   --> Hole: just for constant r
s.Sampling  = [];                 % 'Isolated'/'Oversampling'
s.Option    = 'Dataset';          % 'Dataset'/'NN'/'HO'/ 'Hybrid'
s.nelem     =  20;                %  Mesh refining
s.Print     = false;

% UNIFORM DISTRIBUTION
s.r = ones(5,15)*0.4
% s.r = [0.1,0.2,0.3,0.4,0.5
%          0.1,0.2,0.3,0.4,0.5
%          0.1,0.2,0.3,0.4,0.5];

%NON-UNIFORM DISTRIBUTION
% s.r=[0.1,0.2,0.3
%      0.4,0.5,0.6
%      0.7,0.8,0.8];


%% Multiscale Isolated
s.Training  = 'Multiscale';
s.Sampling  = 'Isolated';
Mult= CoarseTesting_2D(s);
Mult.compute();

%% EIFEM Isolated
s.Training  = 'EIFEM';
s.Sampling  = 'Isolated';
EIFE_IS= CoarseTesting_2D(s);
EIFE_IS.compute();

%% EIFEM Oversampling
s.Training  = 'EIFEM';
s.Sampling  = 'Oversampling';
EIFE_OV= CoarseTesting_2D(s);
EIFE_OV.compute();

%% COMPARE PLOTS

pos= [545   315   914   498];
% RESIDUAL
figure
set(gcf, 'Position', pos) 
plot(Mult.residualCG,'linewidth',2)
hold on
plot(Mult.residualILU,'linewidth',2)
hold on
plot(Mult.residualPCG,'linewidth',2)
hold on
plot(EIFE_IS.residualPCG,'linewidth',2)
hold on
plot(EIFE_OV.residualPCG,'linewidth',2)

set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")
legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-EIFEM(Isolated)-ILU','ILU-EIFEM(Oversampling)-ILU'});


%% 3D case

% case parameters
s.Training  = 'EIFEM';                 % 'EIFEM'/'Multiscale'
s.Sampling  = 'Oversampling';          % 'Isolated'/'Oversampling'
s.Option    = 'Direct';          % % 'Dataset'/'NN'/'Direct

% UNIFORM DISTRIBUTION
s.r = ones(3,7,3)*0.4;
s.nelem     =  10;                %  Mesh refining

s.fileNameEIFEM = 'Sphere_r04_Over.mat';
t_EIFE= CoarseTesting_3D(s);
t_EIFE.compute();


s.fileNameEIFEM = 'Sphere_r04_Mult.mat';
t_Over= CoarseTesting_3D(s);
t_Over.compute();


figure
set(gcf, 'Position', pos) 
plot(t_EIFE.residualCG,'linewidth',2)
hold on
plot(t_EIFE.residualILU,'linewidth',2)
hold on
plot(t_Over.residualPCG,'linewidth',2)
hold on
plot(t_EIFE.residualPCG,'linewidth',2)
set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")

legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-EIFEM(dirac)-ILU'});