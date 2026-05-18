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
s.fileNameEIFEM = [];

% UNIFORM DISTRIBUTION
% s.r = ones(3,10)*0.8;

% NON-UNIFORM DISTRIBUTION
s.r= [ 0.25, 0.40, 0.55, 0.20, 0.45, 0.60, 0.35, 0.25, 0.50, 0.30;
       0.50, 0.20, 0.35, 0.60, 0.25, 0.40, 0.45, 0.55, 0.30, 0.20;
       0.35, 0.55, 0.45, 0.30, 0.50, 0.25, 0.60, 0.40, 0.20, 0.55];

% OUT DATASET
% s.r= [ 0.275, 0.40, 0.575, 0.2250, 0.45, 0.6250, 0.375, 0.275, 0.54, 0.3250;
%        0.530, 0.2750, 0.375, 0.6250, 0.275, 0.430, 0.475, 0.525, 0.3250, 0.2250;
%        0.375, 0.575, 0.425, 0.3250, 0.5250, 0.275, 0.6250, 0.430, 0.230, 0.575];


% MULTISCALE ISOLATED
s.Training  = 'Multiscale';
s.Sampling  = 'Isolated';
Mult= CoarseTesting_2D(s);
Mult.compute();
% 
% % MULTISCALE ISOLATED + NN
% s.Option = 'NN';
% Mult_NN  = CoarseTesting_2D(s);
% Mult_NN.compute();

% EIFEM ISOLATED
% s.Training  = 'EIFEM';
% s.Sampling  = 'Isolated';
% s.Option    = 'Dataset';
% EIFE_IS= CoarseTesting_2D(s);
% EIFE_IS.compute();


% EIFEM OVERSAMPLING
s.Training  = 'EIFEM';  
s.Sampling  = 'Oversampling';
EIFE_OV= CoarseTesting_2D(s);
EIFE_OV.compute();


% COMPARE PLOTS
pos= [545   315   673   498];
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
legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-EIFEM(Isolated)-ILU','ILU-EIFEM(Oversampling)-ILU'});


%% 2D CASE - 1 PARAM CIRCLE ITER VS RADIUS

s.Training  = [];                 % 'EIFEM'/'Multiscale'
s.Inclusion = 'Material';         % 'Hole'/'Material'/'HoleRaul'   --> Hole: just for constant r
s.Sampling  = [];                 % 'Isolated'/'Oversampling'
s.Option    = 'Dataset';          % 'Dataset'/'NN'/'HO'/ 'Hybrid'
s.nelem     =  20;                %  Mesh refining
s.Print     = false;

r= 0.05:0.05:0.9;

iter_ILU=zeros(1,length(r));
iter_Mult=zeros(1,length(r));
iter_NN=zeros(1,length(r));
iter_Iso=zeros(1,length(r));
iter_Over=zeros(1,length(r));

for i=1:length(r)
    s.r = ones(3,10)*r(i);

    % MULTISCALE ISOLATED
    s.Training  = 'Multiscale';
    s.Sampling  = 'Isolated';
    Mult= CoarseTesting_2D(s);
    Mult.compute();


    % MULTISCALE ISOLATED + NN
    s.Option = 'NN';
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

    iter_ILU(i)=size(Mult.residualILU,2);
    iter_Mult(i)=size(Mult.residualPCG,2);
    iter_NN(i)=size(Mult_NN.residualPCG,2);
    iter_Iso(i)=size(EIFE_IS.residualPCG,2);
    iter_Over(i)=size(EIFE_OV.residualPCG,2);
end

% plot iter vs r
figure
plot (r,iter_ILU, 'LineWidth',1.5);
hold on
plot (r,iter_Mult,'LineWidth',1.5);
hold on
plot (r,iter_NN,'LineWidth',1.5);
hold on,
plot(r,iter_Iso,'LineWidth',1.5);
hold on
plot(r,iter_Over,'LineWidth',1.5);

title('Number of iterations vs radius');
xlabel('r');
ylabel('Iterations');
legend({'ILU', 'ILU-Multiscale-ILU','ILU-Multiscale with NN-ILU','ILU-EIFEM(Isolated)-ILU','ILU-EIFEM(Oversampling)-ILU'});
ylim([0 500]);
xlim([0.05 0.9]);

%% DIRECT
s.Training  = [];             % 'EIFEM'/'Multiscale'/'EIFisol'
s.Sampling  = [];           % 'Isolated'/'Oversampling'
s.Inclusion = 'Hole';       % 'Hole'/'Material'/  --> Hole: just for imported meshes or constant geometry
s.Option    = 'Direct';     % 'Dataset'/'NN'/'HO'/ 'Hybrid'/'Direct
s.r         = ones(3,10)*0.8;
s.nelem     = [];            %  Mesh refining
s.fileNameEIFEM = 'r_0.8_test_interior_refinement.mat';
% s.fileNameEIFEM = 'r0.8_RaulOversampling';
% s.fileNameEIFEM = 'r0.8_RaulMultiscale';
testDirect  = CoarseTesting_2D(s);
testDirect.compute();


%% 2D CASE - 2 PARAMS LATTICE

% case parameters
s.Training  = [];                 % 'EIFEM'/'Multiscale'/'EIFisol'
s.Option    = 'Dataset';          % 'Dataset'/'NN'
s.nelem     =  30;                %  Mesh refining

% % UNIFORM DISTRIBUTION
s.tFrame = ones(6,30)*0.25;
s.tCross = ones(6,30)*0.4;
% 

%NON-UNIFORM DISTRIBUTION
% s.tFrame= [0.15, 0.30, 0.45, 0.10, 0.25, 0.40, 0.35, 0.20, 0.15, 0.45;
%            0.40, 0.10, 0.25, 0.35, 0.15, 0.45, 0.20, 0.30, 0.40, 0.10;
%            0.25, 0.45, 0.15, 0.20, 0.40, 0.30, 0.10, 0.35, 0.25, 0.35];
% 
% s.tCross=[0.10, 0.45, 0.60, 0.25, 0.35, 0.50, 0.15, 0.40, 0.20, 0.55;
%           0.50, 0.20, 0.30, 0.60, 0.10, 0.45, 0.55, 0.15, 0.35, 0.25;
%           0.35, 0.55, 0.15, 0.40, 0.60, 0.20, 0.45, 0.30, 0.50, 0.10];

% MULTISCALE ISOLATED
s.Training  = 'Multiscale';
s.Sampling  = 'Isolated';
% Mult= CoarseTesting_2Params(s);
% Mult.compute();

% MULTISCALE + NN
s.Option = 'NN';
Mult_NN  = CoarseTesting_2Params(s);
Mult_NN.compute();


% EIFEM ISOLATED
s.Training  = 'EIFEM';
s.Sampling  = 'Isolated';
s.Option    = 'Dataset';
EIFE_IS= CoarseTesting_2Params(s);
EIFE_IS.compute();


% EIFEM OVERSAMPLING
s.Training  = 'EIFEM';
EIFE_OV= CoarseTesting_2Params(s);
EIFE_OV.compute();


% COMPARE PLOTS
pos= [545   315   673   498];

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
legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-Multiscale with NN-ILU','ILU-EIFEM(Isolated)-ILU','ILU-EIFEM(Oversampling)-ILU'});



%% 3D CASE - SPHERE

pos= [545   315   673   498];
% case parameters
s.Training  = 'Multiscale';       % 'EIFEM'/'Multiscale'/'EIFisol'
s.Inclusion = 'Material';         % 'Hole'/'Material'/  --> Hole: just for imported meshes or constant geometry
s.Option    = 'Dataset';          % % 'Dataset'/'NN'/'Direct
s.Geometry  = 'Sphere';
s.fileNameEIFEM = [];
s.nelem     =  15;                %  Mesh refining

% UNIFORM DISTRIBUTION
s.r = ones(3,8,3)*0.4;

% s.fileNameEIFEM = 'Sphere_r04_Mult.mat';

% NON UNIFORM DISTRIBUTION
% s.r= zeros(2,6,2);
% s.r(:,:,1) = [
%     0.25, 0.45, 0.60, 0.35, 0.20, 0.55;
%     0.50, 0.30, 0.25, 0.65, 0.40, 0.35
% ];
% 
% s.r(:,:,2) = [
%     0.60, 0.20, 0.45, 0.55, 0.30, 0.65;
%     0.35, 0.50, 0.25, 0.40, 0.60, 0.25
% ];


% Multiscale
Sph_Mult= CoarseTesting_3D(s);
Sph_Mult.compute();

s.Option    = 'NN';   
Sph_NN= CoarseTesting_3D(s);
Sph_NN.compute();


% EIFEM Isolated
s.Option    = 'Dataset';  
s.Training  = 'EIFisol';
Sph_Iso= CoarseTesting_3D(s);
Sph_Iso.compute();


% EIFEM Oversampling
s.Training  = 'EIFEM';
Sph_Over= CoarseTesting_3D(s);
Sph_Over.compute();



figure
set(gcf, 'Position', pos) 
plot(Sph_Mult.residualCG,'linewidth',2)
hold on
plot(Sph_Mult.residualILU,'linewidth',2)
hold on
plot(Sph_Mult.residualPCG,'linewidth',2)
hold on
plot(Sph_NN.residualPCG,'linewidth',2)
hold on
plot(Sph_Iso.residualPCG,'linewidth',2)
hold on
plot(Sph_Over.residualPCG,'linewidth',2)
set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")

legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-Multiscale with NN-ILU','ILU-EIFEM(Isolated)-ILU','ILU-EIFEM(Oversampling)-ILU'});
% ylim([10^-9 10^3]);
% xlim([0 4061]);


%% 3D CASE - Cube

pos= [545   315   673   498];
% case parameters
s.Training  = 'Multiscale';       % 'EIFEM'/'Multiscale'/'EIFisol'
s.Inclusion = 'Material';         % 'Hole'/'Material'/  --> Hole: just for imported meshes or constant geometry
s.Option    = 'Dataset';          % % 'Dataset'/'NN'/'Direct
s.Geometry  = 'Cube';
s.fileNameEIFEM = [];
s.nelem     =  20;                %  Mesh refining

% UNIFORM DISTRIBUTION
s.r = ones(8,4,3)*1;

% NON UNIFORM DISTRIBUTION
% s.r= zeros(2,6,2);
% s.r(:,:,1) = [
%     0.6, 1.2, 0.5, 1.7, 0.9, 1.4;
%     1.1, 0.8, 1.5, 0.7, 1.3, 1.0];
% 
% s.r(:,:,2) = [
%    1.6, 0.5, 1.1, 1.4, 0.8, 1.7;
%    0.9, 1.3, 0.6, 1.2, 1.5, 0.7];


% Multiscale
Cub_Mult= CoarseTesting_3D(s); 
Cub_Mult.compute();

s.Option    = 'NN';   
Cub_NN= CoarseTesting_3D(s);
Cub_NN.compute();


% EIFEM Isolated
s.Option    = 'Dataset';  
s.Training  = 'EIFisol';
Cub_Iso= CoarseTesting_3D(s);
Cub_Iso.compute();


% EIFEM Oversampling
s.Training  = 'EIFEM';
Cub_Over= CoarseTesting_3D(s);
Cub_Over.compute();



figure
set(gcf, 'Position', pos) 
plot(Cub_Mult.residualCG,'linewidth',2)
hold on
plot(Cub_Mult.residualILU,'linewidth',2)
hold on
plot(Cub_Mult.residualPCG,'linewidth',2)
hold on
plot(Cub_NN.residualPCG,'linewidth',2)
hold on
plot(Cub_Iso.residualPCG,'linewidth',2)
hold on
plot(Cub_Over.residualPCG,'linewidth',2)
set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")

legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-Multiscale with NN-ILU','ILU-EIFEM(Isolated)-ILU','ILU-EIFEM(Oversampling)-ILU'});
% ylim([10^-9 10^3]);
% xlim([0 4061]);



%% 3D CASE - Lattice

pos= [545   315   673   498];
% case parameters
s.Training  = 'Multiscale';       % 'EIFEM'/'Multiscale'/'EIFisol'
s.Inclusion = 'Material';         % 'Hole'/'Material'/  --> Hole: just for imported meshes or constant geometry
s.Option    = 'Dataset';          % % 'Dataset'/'NN'/'Direct
s.Geometry  = 'Lattice3D';
s.fileNameEIFEM = [];
s.nelem     =  10;                %  Mesh refining

% UNIFORM DISTRIBUTION
s.tFrame = ones(6,2,2)*0.25;
s.tCross = ones(6,2,2)*0.2;

% NON UNIFORM DISTRIBUTION
% s.r= zeros(2,6,2);
% s.r(:,:,1) = [
%     0.6, 1.2, 0.5, 1.7, 0.9, 1.4;
%     1.1, 0.8, 1.5, 0.7, 1.3, 1.0];
% 
% s.r(:,:,2) = [
%    1.6, 0.5, 1.1, 1.4, 0.8, 1.7;
%    0.9, 1.3, 0.6, 1.2, 1.5, 0.7];


% Multiscale
Latt_Mult= CoarseTesting_3D(s); 
Latt_Mult.compute();

s.Option    = 'NN';   
Latt_NN= CoarseTesting_3D(s);
Latt_NN.compute();


% EIFEM Isolated
s.Option    = 'Dataset';  
s.Training  = 'EIFisol';
Latt_Iso= CoarseTesting_3D(s);
Latt_Iso.compute();


% EIFEM Oversampling
s.Training  = 'EIFEM';
Latt_Over= CoarseTesting_3D(s);
Latt_Over.compute();



figure
set(gcf, 'Position', pos) 
plot(Cub_Mult.residualCG,'linewidth',2)
hold on
plot(Cub_Mult.residualILU,'linewidth',2)
hold on
plot(Cub_Mult.residualPCG,'linewidth',2)
hold on
plot(Cub_NN.residualPCG,'linewidth',2)
hold on
plot(Cub_Iso.residualPCG,'linewidth',2)
hold on
plot(Cub_Over.residualPCG,'linewidth',2)
set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")

legend({'CG', 'ILU', 'ILU-Multiscale-ILU','ILU-Multiscale with NN-ILU','ILU-EIFEM(Isolated)-ILU','ILU-EIFEM(Oversampling)-ILU'});
% ylim([10^-9 10^3]);
% xlim([0 4061]);


%% 3D CASE - AIRFOIL


% case parameters
s.Training  = [];             % 'EIFEM'/'Multiscale'/'EIFisol'
s.Sampling  = [];           % 'Isolated'/'Oversampling'
s.Inclusion = 'Hole';       % 'Hole'/'Material'/  --> Hole: just for imported meshes or constant geometry
s.Option    = 'Direct';     % 'Dataset'/'NN'/'Direct
s.r         = [];
s.nelem     = [];            %  Mesh refining
s.Geometry  = [];

% s.fileNameEIFEM = 'Airfoil_Isolated.mat';
% Air_Iso= CoarseTesting_3D(s);
% Air_Iso.compute();

% 
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

pos= [545   315   673   498];
figure
set(gcf, 'Position', pos) 

% plot(Air_Mult.residualCG,'linewidth',2)
% hold on
plot(Air_Mult.residualILU,'linewidth',2)
hold on
plot(Air_Mult.residualPCG,'linewidth',2)
hold on
plot(Air_ref.residualPCG,'linewidth',2)

set(gca, 'YScale', 'log')
xlabel('Iteration')
ylabel('Residual')
title("Residual evolution")

legend({'ILU','ILU-Multiscale-ILU','ILU-EIFEM(Oversampling)-ILU'});
% ylim([0,37340])