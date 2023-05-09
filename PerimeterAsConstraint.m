filename = 'SquareForAniTestsSymmetric';
% filename = 'CantileverVertical';
% filename = 'CantileverVerticalSymmetric';
% filename = 'CantileverBeam_Triangle_Linear';

ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';

initial_case = 'circleInclusion'; % squareInclusion

cost = {'chomog_alphabeta','anisotropicPerimeterInterior2D'}; % anisotropicPerimeter2D
weights = [1,0]; % 0.21, 0.30

constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};

optimizerUnconstrained = 'SLERP'; 
optimizer = 'NullSpace';
incrementFactor = 2;
designVariable = 'LevelSet';
filterType = 'P1';
fracRadius = 0.4;
micro.alpha =[1 1 0]';
micro.beta =[1 1 0]';

nsteps = 1; % 200
Vfrac_final = 0.6;
optimality_final =1e-3;
constr_final =1e-3;

% Vfrac_initial = 1-pi*(fracRadius/2)^2;
Vfrac_initial = 0.6;
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
plotting = false;
printing = true;
printing_physics = false;
monitoring = true;
monitoring_interval = 1;
maxiter = 400; % 2100   300*nsteps