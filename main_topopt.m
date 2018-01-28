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
settings=struct;
settings.filename='TOPOPT_TEST';
%settings.filename='topopt_quad';
%settings.filename='GrippingNew';

settings.plotting=true;
settings.printing=false;
settings.maxiter = 5000;

settings.method='SIMPALL';
%settings.method='SIMP_P3';
% settings.method='SIMP_Adaptative';

settings.material='ISOTROPIC';
settings.initial_case='full';

%settings.ptype='Compliance_st_Volume';
settings.ptype='ComplianceLamPerimeter_st_Volume';
%settings.ptype='Compliance_st_VolumePerimeter';

if settings.filename=='GrippingNew'
    settings.ptype='Gripping';
end

settings.optimizer='SLERP';
%settings.optimizer='PROJECTED GRADIENT';
%settings.optimizer='MMA';
%settings.optimizer='IPOPT';

settings.filter='P1';
%settings.filter='PDE';
settings.TOL.rho_plus=1;
settings.TOL.rho_minus=0;
settings.TOL.E_plus=1;
settings.TOL.E_minus=1e-3;
settings.TOL.nu_plus=1/3;
settings.TOL.nu_minus=1/3;

settings.target_parameters.Vfrac=0.3;
settings.target_parameters.optimality_tol=1e-3;
settings.target_parameters.constr_tol=1e-3;
settings.target_parameters.Perimeter_target=5;
settings.perimeter.optimizer=settings.optimizer;
settings.perimeter.lambda=0.1;

settings.nsteps=1;
settings.Vfrac_final=settings.target_parameters.Vfrac;
settings.optimality_final=settings.target_parameters.optimality_tol;
settings.constr_final=settings.target_parameters.constr_tol;
settings.Vfrac_initial=1;
settings.optimality_initial=1e-3;
settings.constr_initial=1e-3;
%% main
tic
test=TopOpt_Problem.create(settings);
test.preProcess;
test.computeVariables;
toc
%% Video creation
if settings.printing
gidPath = 'C:\Program Files\GiD\GiD 13.0.2';% 'C:\Program Files\GiD\GiD 13.0.3';
files_name = [];
files_folder = fullfile(pwd,'Output');
iterations = 0:test.optimizer.niter;
video_name=strcat('2Video_',settings.ptype,'_',settings.optimizer,'_',settings.method,'_',int2str(settings.nsteps),'_0dot',int2str(10*settings.Vfrac_final),'.gif');
My_VideoMaker = VideoMaker_TopOpt.Create(settings.optimizer);
My_VideoMaker.Set_up_make_video(gidPath,files_name,files_folder,iterations)
% 
output_video_name_design_variable = fullfile(pwd,video_name);
My_VideoMaker.Make_video_design_variable(output_video_name_design_variable)

% %
% output_video_name_design_variable_reg = fullfile(pwd,'DesignVariable_Reg_Video');
% My_VideoMaker.Make_video_design_variable_reg(output_video_name_design_variable_reg)
% 
% output_video_name_design_variable_reg = fullfile(pwd,'DesignVariable_Reg_Video');
% My_VideoMaker.Make_video_design_variable_reg(output_video_name_design_variable_reg)
% 
% output_video_name_stress = fullfile(pwd,'Stress_Video');
% My_VideoMaker.Make_video_stress(output_video_name_stress)
end
