% filename = 'SquareForAniTests';
filename = 'CantileverVertical';
% filename = 'CantileverBeam_Triangle_Linear';

ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';

initial_case = 'full'; % squareInclusion

cost = {'compliance'}; % anisotropicPerimeterInterior2D
weights = [1];

constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};

optimizerUnconstrained = 'SLERP'; 
optimizer = 'NullSpace';
incrementFactor = 2;
designVariable = 'LevelSet';
filterType = 'P1';
widthSquare = sqrt(0.15);

nsteps = 1;
% Vfrac_final = 0.85;
Vfrac_final = 0.4;
optimality_final =1e-3;
constr_final =1e-3;

% Vfrac_initial = 0.85;
Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;
Perimeter_target = 5;
perimeterTarget = 5;

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
maxiter = 500; % 3000    200      50