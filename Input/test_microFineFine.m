filename='RVE_Square_Triangle_FineFine';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'squareInclusion';
cost={'anisotropicPerimeterInterior2D'};
weights=[1];
constraint = {'volumeConstraint'};
optimizerUnconstrained = 'PROJECTED GRADIENT';
optimizer = 'DualNestedInPrimal';
incrementFactor = 1.5;
designVariable = 'Density';
filterType = 'P1';

widthSquare = 0.5;

nsteps = 15;
Vfrac_final = 0.75;
Perimeter_target=5;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.75;
optimality_initial = 1e-3;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

%Micro
epsilon_isotropy_initial=1e-1;
epsilon_isotropy_final = 1e-3;
micro.alpha =[1 1 0]';
micro.beta =[1 1 0]';

% For all tests
plotting = true;
printing = false;
monitoring = false;
monitoring_interval = 1;
maxiter = 900;
