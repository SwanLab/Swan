clc
clear variables;close all;
addpath(genpath('.\FEM\'));
addpath(genpath('.\Topology Optimization\'));
%% test
run('test.m');
clear variables;
%% settings
settings=struct;
settings.filename='TOPOPT_TEST';

% settings.method='SIMP_P3';
settings.method='SIMP_Adaptative';
% settings.method='SIMPALL';

settings.material='ISOTROPIC';
settings.ptype='Compliance_st_Volume';
settings.initial_case='full';

settings.optimizer='SLERP';
%settings.optimizer='PROJECTED GRADIENT';
%settings.optimizer='MMA';

settings.filter='P1';
settings.TOL.rho_plus=1;
settings.TOL.rho_minus=0;
settings.TOL.E_plus=1;
settings.TOL.E_minus=1e-3;
settings.TOL.nu_plus=1/3;
settings.TOL.nu_minus=1/3;
settings.epsilon_scalar_product_P1=0.03;
settings.volume.Vfrac=0.3;
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