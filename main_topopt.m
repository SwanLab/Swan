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

%settings.filename='TOPOPT_TEST';  %MACRO
%settings.filename='RVE_Square_Triangle';
%settings.filename='RVE_Square_Triangle_Fine';
settings.filename='topopt_quad';
%settings.filename='GrippingNew';

settings.plotting=true;
settings.printing=false;
settings.maxiter = 5000;


settings.method='SIMPALL';
%settings.method='SIMP_P3';
% settings.method='SIMP_Adaptative';

settings.material='ISOTROPIC';
settings.initial_case='full';
%settings.initial_case='circle';
% settings.initial_case='horizontal';
% settings.initial_case='square';
% settings.initial_case='feasible';
% settings.initial_case='rand';


settings.ptype='Compliance_st_Volume';
%settings.ptype='ComplianceLamPerimeter_st_Volume';
% settings.ptype='Compliance_st_VolumePerimeter';
% settings.ptype='Chomog_alphabeta_st_Volume';
% settings.ptype='Chomog_fraction_st_Volume';
%settings.ptype='ChomogLamPerimeter_alphabeta_st_Volume';
%settings.ptype='ChomogLamPerimeter_fraction_st_Volume';

%if settings.filename=='GrippingNew'
%    settings.ptype='Gripping';
%end

%settings.optimizer='SLERP';
settings.optimizer='PROJECTED GRADIENT';
%settings.optimizer='MMA';
%settings.optimizer='IPOPT';


settings.filter='P1'; %PDE

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
settings.perimeter.lambda=0.1;%%%%%%%%%%%%

settings.nsteps=1;
settings.Vfrac_final=settings.target_parameters.Vfrac;
settings.optimality_final=settings.target_parameters.optimality_tol;
settings.constr_final=settings.target_parameters.constr_tol;
settings.Vfrac_initial=1;

settings.optimality_initial=1e-1;
settings.constr_initial=1e-1;


settings.micro.alpha =[1 0 0]';
settings.micro.beta =[1 0 0]';


settings.optimality_initial=1e-3;
settings.constr_initial=1e-3;

%% main

tic
test=TopOpt_Problem.create(settings);
% test.checkDerivative;
% toc
test.preProcess;
test.computeVariables;
toc

test.postProcess;


