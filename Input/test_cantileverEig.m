filename = 'CantileverBeam_Triangle_Linear';%'Cantilever';%'CantileverBeam_Triangle_Linear';
ptype = 'MACRO';
method = 'SIMP_P3';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'firstEigTopCost'};
weights = [1];
constraint = {'ModalfirstEigConstraint','volumeConstraint'};
constraint_case = {'INEQUALITY','INEQUALITY'};
optimizerUnconstrained = 'PROJECTED GRADIENT';%'SLERP';%'PROJECTED GRADIENT'; 
optimizer = 'MMA'; %'DualNestedInPrimal';'DualNestedInPrimal';%'AlternatingPrimalDual';%'AlternatingPrimalDual';
incrementFactor = 1;
designVariable = 'DensityEigModes';%'LevelSet';%'Density';
filterType = 'P1';

nsteps = 1;
Vfrac_final = 0.4;
optimality_final =1e-5;
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
maxiter = 1000;