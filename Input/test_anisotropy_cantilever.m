% filename = 'jaCantilever';
% filename = 'Bridge_UltraFine';
filename = 'ArchUltraFine';
% filename = 'MicroUltraFine';

%Micro
% epsilon_isotropy_initial=1e-1;
% epsilon_isotropy_final = 1e-3;
% micro.alpha =[1 1 0]';
% micro.beta =[1 1 0]';

ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
% constraint_case = 'EQUALITY';
optimizerUnconstrained = 'SLERP';
optimizer = 'DualNestedInPrimal';
incrementFactor = 1.5; % Recommended: 1.5; 2.0
designVariable = 'LevelSet';
filterType = 'P1';

nsteps = 4; % Recommended slope: 5%/step is OK
Vfrac_final = 0.15;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;
Perimeter_target = 5;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

% For all tests
plotting = true;
printing = false;
printing_physics = false;
monitoring = true;
monitoring_interval = 1;
maxiter = 2000;