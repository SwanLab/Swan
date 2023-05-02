filename='test3d_micro_cube';
% filename='holeinclusion3d';
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'sphereInclusion';
% cost={'chomog_alphabeta','perimeterConstraint'};
cost={'chomog_alphabeta'};
weights=[1];
constraint = {'volumeConstraint'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizerUnconstrained = 'SLERP'; 
% optimizer = 'DualNestedInPrimal';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizerUnconstrained = 'PROJECTED GRADIENT'; 
% optimizer = 'NullSpace';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimizer = 'MMA';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
incrementFactor = 1;
designVariable = 'Density';
filterType = 'P1';
fracRadius = 0.50; %%%%%%%%%%% 0.75

nsteps = 1;
Vfrac_final = 0.5;
Perimeter_target=1;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.75; %%%%%%%%%%% 1
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
micro.alpha =[1 1 0 0 0 0]';
micro.beta =[1 1 0 0 0 0]';

% For all tests
plotting = false;
printing = true;
monitoring = true;
maxiter = 150;
monitoring_interval = 1;
