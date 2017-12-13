clc
clear variables;close all;
 
addpath(genpath(fullfile('.','FEM')));
addpath(genpath(fullfile('.','Topology Optimization')));
%% test
run('test.m');
clear variables;
%% settings
settings=struct;
settings.filename='TOPOPT_TEST';

% settings.method='SIMP_P3';
% settings.method='SIMP_Adaptative';
settings.method='SIMPALL';

settings.material='ISOTROPIC';
settings.ptype='Compliance_st_Volume';
settings.initial_case='full';

settings.optimizer='SLERP';
%settings.optimizer='PROJECTED GRADIENT';
% settings.optimizer='MMA';

settings.filter='P1';
settings.TOL.rho_plus=1;
settings.TOL.rho_minus=0;
settings.TOL.E_plus=1;
settings.TOL.E_minus=1e-3;
settings.TOL.nu_plus=1/3;
settings.TOL.nu_minus=1/3;
settings.epsilon_scalar_product_P1=0.03;
settings.volume.Vfrac=0.4;
%% main
tic
switch settings.ptype
    case 'Compliance_st_Volume'
        settings.nconstr=1;
        test=TopOpt_Problem_Compliance_st_Volume(settings);
    otherwise
        disp('Problem not added')
end
test.preProcess;
test.computeVariables;
toc
%% Video creation
post = Postprocess;
gidPath = '/opt/GiDx64/13.0.2/';
files_name = test.settings.filename;
files_folder = fullfile(pwd,'Output','TOPOPT_RESULTS');
iterations = 1:test.optimizer.niter;

output_video_name = fullfile(pwd,'CharacteristicVideo');
post.Print_make_video_characteristic_function(gidPath,files_name,files_folder,iterations,output_video_name)


output_video_name = fullfile(pwd,'StressVideo');
post.Print_make_video_stress(gidPath,files_name,files_folder,iterations,output_video_name)