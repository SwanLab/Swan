filename = 'CantileverVertical';

ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';

initial_case = 'full';

cost = {'compliance','anisotropicPerimeter2D'};
weights = [1,1e-5];

constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};

optimizerUnconstrained = 'SLERP'; 
optimizer = 'AlternatingPrimalDual';
incrementFactor = 2;
designVariable = 'LevelSet';
filterType = 'P1';

nsteps = 11;
Vfrac_final = 0.4;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.4; % 0.4 !!!!!
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
monitoring_interval = 5;
maxiter = 2200;