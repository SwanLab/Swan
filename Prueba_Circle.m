% Example
clear variables;close all;
 
addpath(genpath('./FEM'));
addpath(genpath('./Topology Optimization'));

settings.filename='TOPOPT_TEST_MICRO';
settings.ptype='Chomog_alphabeta';
settings.method='SIMPALL';
settings.material='ISOTROPIC';
settings.optimizer='SLERP';
settings.initial_case='circle';
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
settings.maxiter = 5;

settings.micro.alpha =sqrt(2)/2*[1 1 0]';
settings.micro.beta =sqrt(2)/2*[1 1 0]';

%% main

Micro_Square_Triangle = Physical_Problem_Micro('RVE_Square_Triangle');
Micro_Square_Triangle.preProcess;

test=TopOpt_Problem.create(settings);
test.preProcess;

rho_nodal = test.x;
rho_elem = test.filter.getP0fromP1(rho_nodal);
matprop = test.interpolation.computeMatProp(rho_elem);


Micro_Square_Triangle.setMatProps(matprop)
Micro_Square_Triangle.computeChomog;

Micro_Square_Triangle.variables.Chomog
fprintf('Ok\n');

