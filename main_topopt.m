clc
clear variables;close all;

addpath(genpath('./FEM'));
addpath(genpath('./Topology Optimization'));
% addpath(genpath(fullfile('.','FEM')));
% addpath(genpath(fullfile('.','Topology Optimization')));
%% test
%run('test.m');
clear variables;
%% settings
% settings.filename='TOPOPT_TEST';  %MACRO
settings.filename='RVE_Square_Triangle';
% settings.filename='RVE_Square_Triangle_Fine';



% settings.method='SIMP_P3';
% settings.method='SIMP_Adaptative';
settings.method='SIMPALL';

settings.material='ISOTROPIC';
% settings.initial_case='full';
settings.initial_case='circle';
% settings.initial_case='horizontal';
% settings.initial_case='square';
% settings.initial_case='feasible';
% settings.initial_case='rand';



% settings.ptype='Compliance_st_Volume';
% settings.ptype='ComplianceLamPerimeter_st_Volume';
% settings.ptype='Compliance_st_VolumePerimeter';
% settings.ptype='Chomog_alphabeta_st_Volume';
% settings.ptype='Chomog_fraction_st_Volume';
settings.ptype='ChomogLamPerimeter_alphabeta_st_Volume';
% settings.ptype='ChomogLamPerimeter_fraction_st_Volume';



settings.optimizer='SLERP';
% settings.optimizer='PROJECTED GRADIENT';
% settings.optimizer='MMA';
% settings.optimizer='IPOPT';

settings.filter='P1'; %PDE
settings.TOL.rho_plus=1;
settings.TOL.rho_minus=0;
settings.TOL.E_plus=1;
settings.TOL.E_minus=1e-3;
settings.TOL.nu_plus=1/3;
settings.TOL.nu_minus=1/3;

settings.target_parameters.Vfrac=0.5;
settings.target_parameters.optimality_tol=1e-3;
settings.target_parameters.constr_tol=1e-3;
settings.target_parameters.Perimeter_target=3.5;
settings.perimeter.optimizer=settings.optimizer;
settings.perimeter.lambda=0.1;%%%%%%%%%%%%

settings.nsteps=1;
settings.Vfrac_final=settings.target_parameters.Vfrac;
settings.optimality_final=settings.target_parameters.optimality_tol;
settings.constr_final=settings.target_parameters.constr_tol;
settings.Vfrac_initial=1;
settings.optimality_initial=1e-1;
settings.constr_initial=1e-1;
settings.maxiter = 200;

settings.micro.alpha =[1 1 0]';
settings.micro.beta =[1 1 0]';

%% main
tic
test=TopOpt_Problem.create(settings);
test.preProcess;
test.computeVariables;
toc
%% Video creation
gidPath = 'C:\Program Files\GiD\GiD 13.0.3';
files_name = test.settings.filename;
files_folder = fullfile(pwd,'Output');
iterations = 1:test.optimizer.niter;

My_VideoMaker = VideoMaker_TopOpt.Create(settings.optimizer);
My_VideoMaker.Set_up_make_video(gidPath,files_name,files_folder,iterations)

% output_video_name_design_variable_reg = fullfile(pwd,'DesignVariable_Reg_Video');
% My_VideoMaker.Make_video_design_variable_reg(output_video_name_design_variable_reg)

output_video_name_design_variable = fullfile(pwd,'DesignVariable_Video');
My_VideoMaker.Make_video_design_variable(output_video_name_design_variable)

% output_video_name_design_variable_reg = fullfile(pwd,'DesignVariable_Reg_Video');
% My_VideoMaker.Make_video_design_variable_reg(output_video_name_design_variable_reg)

% output_video_name_stress = fullfile(pwd,'Stress_Video');
% My_VideoMaker.Make_video_stress(output_video_name_stress)

