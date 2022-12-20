filename = 'anisoCantilever';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
cost = {'compliance','anisotropicPerimeterInterior2D'};
weights = [1,0.05];
constraint = {'volumeConstraint'};
optimizerUnconstrained = 'SLERP';
optimizer = 'DualNestedInPrimal';
incrementFactor = 1.5;
designVariable = 'LevelSet';
filterType = 'P1';
anisoScaleAngle = 60;
anisoOverhangAngle = 90;

nsteps = 1;
Vfrac_final = 0.6;
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
plotting = false;
printing = false;
printing_physics = false;
monitoring = false;
maxiter = 1;