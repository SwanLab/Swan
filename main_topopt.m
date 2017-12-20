clc
clear variables;close all;
addpath(genpath('.\FEM\'));
addpath(genpath('.\Topology Optimization\'));
%% test
%run('test.m')
clear variables
%% settings
settings=struct;
settings.filename='TOPOPT_TEST';
settings.method='SIMPALL';
settings.material='ISOTROPIC';
settings.initial_case='full';

settings.ptype='Compliance_st_Volume';
%settings.ptype='ComplianceLamPerimeter_st_Volume';
%settings.ptype='Compliance_st_VolumePerimeter';

settings.optimizer='SLERP';
%settings.optimizer='PROJECTED GRADIENT';
%settings.optimizer='MMA';
%settings.optimizer='IPOPT';

settings.filter='P1';
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
settings.perimeter.lambda=0.1;

settings.nsteps=1;
settings.Vfrac_final=settings.target_parameters.Vfrac;
settings.optimality_final=settings.target_parameters.optimality_tol;
settings.constr_final=settings.target_parameters.constr_tol;
settings.Vfrac_initial=1;
settings.optimality_initial=1e-1;
settings.constr_initial=1e-1;
%% main
tic
test=TopOpt_Problem.create(settings);
test.preProcess;
test.computeVariables;

toc