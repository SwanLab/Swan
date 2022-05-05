% filename='RVE_Square_Triangle_FineFine';
filename='test2d_micro';
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
cost={'chomog_alphabeta'};
weights=[1];
constraint = {'volumeConstraint'};
optimizerUnconstrained = 'SLERP';
optimizer = 'DualNestedInPrimal';
incrementFactor = 1.5;
designVariable = 'LevelSet';
filterType = 'P1';
            
% levelFibers = 3;
% volumeFibers = 0.5;

nsteps = 1;
Vfrac_final = 0.5;
Perimeter_target=1;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
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
micro.alpha =[1 1 1]';
micro.beta =[1 1 1]';

% For all tests
plotting = true;
printing = false;
monitoring = true;
monitoring_interval = 1;
maxiter = 2000;
