filename = 'CantileverBeam_Triangle_Linear'; %'Cantilever';%'CantileverBeam_Triangle_Linear';
ptype = 'MACRO';
method = 'SIMP_P3';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'ShFunc_firstEigTop'};
weights = [1];
% Create all shape functionals 
constraint = {'firstFreq','undampedFreeVibrationProblem','orthonormalize','volumeConstraint','minDensity'}; % Change names
optimizerUnconstrained = 'PROJECTED GRADIENT'; % 'PROJECTED GRADIENT'; 
optimizer = 'MMA'; %'NullSpace';%'DualNestedInPrimal';'DualNestedInPrimal';%'AlternatingPrimalDual';%'AlternatingPrimalDual';
incrementFactor = 1;
designVariable = 'Density';%'Density';
filterType = 'P1';

nsteps = 3;
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
printing = false;
printing_physics = false;
monitoring = true;
maxiter = 1000;