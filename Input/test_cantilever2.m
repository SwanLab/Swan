filename = 'MBB_Dapogny_fine';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';%'holes';
cost = {'compliance','anisotropicPerimeter2D'};
weights = [1,0.4];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
optimizerUnconstrained = 'SLERP';%'SLERP';%'PROJECTED GRADIENT';%'PROJECTED GRADIENT'; 
optimizer = 'NullSpace';%'DualNestedInPrimal';'DualNestedInPrimal';%'AlternatingPrimalDual';%'AlternatingPrimalDual';
incrementFactor = 1.2;
designVariable = 'LevelSet';%'Density';'LevelSet'
filterType = 'P1';
nsteps = 200;
Vfrac_final = 0.4;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.4;
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
plotting = false;
printing = false;
printing_physics = false;
monitoring = true;
monitoring_interval = 1;
maxiter = nsteps*3;