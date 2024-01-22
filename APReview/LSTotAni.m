cost = {'nonadjoint_compliance','anisotropicPerimeter2D'};
weights = [1,0.05];
optimizerUnconstrained = 'SLERP';
designVariable = 'LevelSet';
nsteps = 100;
maxiter = 400;






















filename = 'grippingTrialFine';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
fracRadius = 0.05;


constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
optimizer = 'NullSpace';
incrementFactor = 1.2;
filterType = 'P1';
Vfrac_final = 0.25;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.25;
optimality_initial = 1e-3;
constr_initial = 1e-3;
Perimeter_target = 5;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;
plotting = false;
printing = false;
printing_physics = false;
monitoring = true;
monitoring_interval = 1;