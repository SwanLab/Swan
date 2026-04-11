%% PRECONDITIONING SCRIPT

%This script is intented to execute the coarseTesting_Abril for varius
%cases and create plots that summarise all the results in a better way

%% 2D CASE - 1 PARAM CIRCLE

% case parameters
s.Training  = [];                 % 'EIFEM'/'Multiscale'
s.Inclusion = 'Material';         % 'Hole'/'Material'/'HoleRaul'   --> Hole: just for constant r
s.Sampling  = [];                 % 'Isolated'/'Oversampling'
s.Option    = 'Dataset';          % 'Dataset'/'NN'/'HO'/ 'Hybrid'
s.nelem     =  20;                %  Mesh refining
s.Print     = false;

% UNIFORM DISTRIBUTION
% s.r = ones(5,15)*0.4
% s.r = [0.1,0.2,0.3,0.4,0.5
%          0.1,0.2,0.3,0.4,0.5
%          0.1,0.2,0.3,0.4,0.5];

%NON-UNIFORM DISTRIBUTION
s.r= [ 0.25, 0.40, 0.55, 0.20, 0.45, 0.60, 0.35, 0.25, 0.50, 0.30;
       0.50, 0.20, 0.35, 0.60, 0.25, 0.40, 0.45, 0.55, 0.30, 0.20;
       0.35, 0.55, 0.45, 0.30, 0.50, 0.25, 0.60, 0.40, 0.20, 0.55];

% MULTISCALE ISOLATED
s.Training  = 'Multiscale';
s.Sampling  = 'Isolated';
Mult= CoarseTesting_2D(s);
Mult.compute();

% MULTISCALE ISOLATED + NN
s.Option = 'NN'
Mult_NN  = CoarseTesting_2D(s);
Mult_NN.compute();

% EIFEM ISOLATED
s.Training  = 'EIFEM';
s.Sampling  = 'Isolated';
s.Option    = 'Dataset';
EIFE_IS= CoarseTesting_2D(s);
EIFE_IS.compute();


% EIFEM OVERSAMPLING
s.Training  = 'EIFEM';
s.Sampling  = 'Oversampling';
EIFE_OV= CoarseTesting_2D(s);
EIFE_OV.compute();


% COMPARE PLOTS
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
plot(Mult_NN.residualPCG,'linewidth',2)
hold on
plot(EIFE_IS.residualPCG,'linewidth',2)
hold on
plot(EIFE_OV.residualPCG,'linewidth',2)

set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")
legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-Multiscale+NN-ILU','ILU-EIFEM(Isolated)-ILU','ILU-EIFEM(Oversampling)-ILU'});



%% 2D CASE - 2 PARAMS LATTICE

% case parameters
s.Training  = [];                 % 'EIFEM'/'Multiscale'/'EIFisol'
s.Option    = 'Dataset';          % 'Dataset'/'NN'
s.nelem     =  30;                %  Mesh refining

% UNIFORM DISTRIBUTION
s.tFrame = ones(3,10)*0.4;
s.tCross = ones(3,10)*0.2;


%NON-UNIFORM DISTRIBUTION
% s.tFrame= [ 0.25, 0.40, 0.55, 0.20, 0.45, 0.60, 0.35, 0.25, 0.50, 0.30;
%        0.50, 0.20, 0.35, 0.60, 0.25, 0.40, 0.45, 0.55, 0.30, 0.20;
%        0.35, 0.55, 0.45, 0.30, 0.50, 0.25, 0.60, 0.40, 0.20, 0.55];

% MULTISCALE ISOLATED
s.Training  = 'Multiscale';
s.Sampling  = 'Isolated';
Mult= CoarseTesting_2Params(s);
Mult.compute();

% MULTISCALE ISOLATED + NN
% s.Option = 'NN'
% Mult_NN  = CoarseTesting_2Params(s);
% Mult_NN.compute();


% EIFEM OVERSAMPLING
s.Training  = 'EIFEM';
EIFE_OV= CoarseTesting_2Params(s);
EIFE_OV.compute();


% COMPARE PLOTS
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
% plot(Mult_NN.residualPCG,'linewidth',2)
% hold on
% plot(EIFE_IS.residualPCG,'linewidth',2)
% hold on
plot(EIFE_OV.residualPCG,'linewidth',2)

set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")
legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-EIFEM(Oversampling)-ILU'});



%% 3D CASE - SPHERE

pos= [545   315   914   498];
% case parameters
s.Training  = 'EIFEM';                 % 'EIFEM'/'Multiscale'/'EIFisol'
s.Sampling  = 'Oversampling';          % 'Isolated'/'Oversampling'
s.Inclusion = 'Material';         % 'Hole'/'Material'/  --> Hole: just for imported meshes or constant geometry
s.Option    = 'Dataset';          % % 'Dataset'/'NN'/'Direct

% UNIFORM DISTRIBUTION
s.r = ones(3,7,3)*0.4;
s.nelem     =  10;                %  Mesh refining

s.fileNameEIFEM = 'Sphere_r04_Over.mat';
Sph_Iso= CoarseTesting_3D(s);
Sph_Iso.compute();


s.fileNameEIFEM = 'Sphere_r04_Mult.mat';
Sph_Mult= CoarseTesting_3D(s);
Sph_Mult.compute();


figure
set(gcf, 'Position', pos) 
plot(Sph_Iso.residualCG,'linewidth',2)
hold on
plot(Sph_Iso.residualILU,'linewidth',2)
hold on
plot(Sph_Mult.residualPCG,'linewidth',2)
hold on
plot(Sph_Iso.residualPCG,'linewidth',2)
set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")

legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-EIFEM(dirac)-ILU'});


%% 3D CASE - AIRFOIL


% case parameters
s.Training  = [];             % 'EIFEM'/'Multiscale'/'EIFisol'
s.Sampling  = [];           % 'Isolated'/'Oversampling'
s.Inclusion = 'Hole';       % 'Hole'/'Material'/  --> Hole: just for imported meshes or constant geometry
s.Option    = 'Direct';     % 'Dataset'/'NN'/'Direct
s.r         = [];
s.nelem     =[];            %  Mesh refining

% s.fileNameEIFEM = 'Airfoil_Isolated.mat';
% Air_Iso= CoarseTesting_3D(s);
% Air_Iso.compute();


s.fileNameEIFEM = 'Airfoil_Multiscale.mat';
Air_Mult= CoarseTesting_3D(s);
Air_Mult.compute();

% s.fileNameEIFEM = 'Airfoil_Oversampling.mat';
% Air_Over= CoarseTesting_3D(s);
% Air_Over.compute();

% s.fileNameEIFEM = 'Airfoil_Over2.mat';
% Air_Over2= CoarseTesting_3D(s);
% Air_Over2.compute();

%%
s.fileNameEIFEM = 'DEF_Q8_wing_1.mat';
Air_ref= CoarseTesting_3D(s);
Air_ref.compute();

pos= [545   315   914   498];
figure
set(gcf, 'Position', pos) 
% plot(Air_Iso.residualCG,'linewidth',2)
% hold on
% plot(Air_Iso.residualILU,'linewidth',2)
% hold on
plot(Air_Mult.residualPCG,'linewidth',2)
hold on
plot(Air_Iso.residualPCG,'linewidth',2)
hold on
plot(Air_Over.residualPCG,'linewidth',2)
hold on
plot(Air_ref.residualPCG,'linewidth',2)

set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")

legend({'Multiscale','EIF-Isolated','EIF-Oversampling','Joaquin'});