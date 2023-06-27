filename = 'MBB_Dapogny_fine';
% filename = 'SquareForAniTests';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';%'holes';
widthSquare = 0.4;
cost = {'compliance','anisotropicPerimeter2DWithPNorm'};
weights = [1,0.75];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
optimizerUnconstrained = 'PROJECTED GRADIENT';%'SLERP';%'PROJECTED GRADIENT';%'PROJECTED GRADIENT'; 
optimizer = 'NullSpace';%'DualNestedInPrimal';'DualNestedInPrimal';%'AlternatingPrimalDual';%'AlternatingPrimalDual';
incrementFactor = 1.2;
designVariable = 'Density';%'Density';'LevelSet'
filterType = 'P1';
nsteps = 150;
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
plotting = true;
printing = false;
printing_physics = false;
monitoring = true;
monitoring_interval = 1;
maxiter = 600;