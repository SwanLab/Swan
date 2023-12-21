%filename='test3d_micro_cube_hexa';
%filename='aba_hexa_mesh_002_matlab';
filename='test3d_micro_cube_hexa_v2';
% filename='test3d_micro_cube_v2';
% filename='test3d_micro_cube';
% filename='holeinclusion3d';
% filename = 'test2d_micro';
%filename = 'aba_unitform_hexa_mesh_0025_MICRO';
%filename = 'aba_unitform_hexa_mesh_003_precise_MICRO';
%filename = 'aba_unitform_hexa_mesh_0033_MICRO';
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
%initial_case = 'sphereInclusion';
initial_case = 'holes';
% cost={'chomog_alphabeta','perimeterConstraint'};
cost={'chomog_alphabeta'};
weights=[1];
constraint = {'volumeConstraint'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizerUnconstrained = 'SLERP'; % level set
% optimizerUnconstrained = 'PROJECTED GRADIENT'; % density
% optimizer = 'DualNestedInPrimal';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimizerUnconstrained = 'SLERP'; 
%optimizerUnconstrained = 'PROJECTED GRADIENT';
optimizer = 'NullSpace';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizer = 'MMA';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
incrementFactor = 1; %%%%%%%%%%% 1.2
designVariable = 'LevelSet';
%designVariable = 'LevelSet';
filterType = 'P1'; % P1 / PDE
fracRadius = 0.4; %%%%%%%%%%% 0.75

kfrac = 1.1;

nsteps = 1; % 50
Vfrac_final = 0.45; %%%%% 0.7
Perimeter_target=1;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1-4/3*pi*(fracRadius/2)^3; %%%%%%%%%%% 1 (primer problema al q intenta convergir)
optimality_initial = 1e-3;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

%Micro
epsilon_isotropy_initial=1e-1;
epsilon_isotropy_final = 1e-3;
%micro.alpha =[1 0 0 0 0 0]';
micro.alpha =[1 0.125 0.125 0.3 0.125 0.3]';
%micro.beta =[1 0 0 0 0 0]';
micro.beta =[1 0.125 0.125 0.3 0.125 0.3]';
% micro.alpha =[0 0 1]';
% micro.beta =[0 0 1]';

% For all tests
plotting = false;
printing = true; %for paraview
monitoring = true;
maxiter = 145;
monitoring_interval = 1;