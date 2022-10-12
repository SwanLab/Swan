filename = 'bridge_10_2';
ptype = 'MACRO';
method = 'SIMP_P3';
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'volume'};
weights = [1];
constraint = {'complianceConstraintC1'};
constraint_case = {'INEQUALITY'};
optimizerUnconstrained = 'PROJECTED GRADIENT';%'PROJECTED GRADIENT'; 
optimizer = 'AlternatingPrimalDual';%'DualNestedInPrimal';'DualNestedInPrimal';%'AlternatingPrimalDual';%'AlternatingPrimalDual';
incrementFactor = 1;
designVariable = 'Density';%'Density';
filterType = 'P1';

nsteps = 1;
Vfrac_final = 0.4;
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
printing = true;
printing_physics = false;
monitoring = true;
monitoring_interval = 1;
maxiter = 1000;
