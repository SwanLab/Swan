filename = 'archExample_80x40';%'cantileverExample_80x40';%'bridgeExample_200x10_case1';%%'archExample_80x40';%'bridgeExample_120x20';%'cantileverExample_120x60';%'cantileverExample_60x40';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'volume'};
weights = 1;
constraint = {'complianceConstraint'};%{'volumeConstraint'};%{'complianceConstraint'};%
constraint_case = {'EQUALITY'};
optimizerUnconstrained = 'PROJECTED GRADIENT'; 
optimizer = 'IPM';%'MMA';%'IPM';%'AlternatingPrimalDual';%'DualNestedInPrimal';%'%'NullSpace';
incrementFactor = 1.05; 
designVariable = 'Density';
filterType = 'P1';
line_search_initiator = 'INCREASING LAST STEP';


nsteps = 1;
Vfrac_final = 0.4;
optimality_final =1e-3;
constr_final = 1e-5;

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

plotting = true;
printing = false;
printing_physics = false;
monitoring = true;
monitoring_interval = 1;
maxiter = 6000;