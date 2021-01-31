filename = 'CantileverBeam_Triangle_Linear_Fine';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
constraint_case = 'EQUALITY';
optimizer = 'DualNestedInPrimal';
optimizerUnconstrained = 'PROJECTED GRADIENT';
incrementFactor = 1;
rate = 0.95;
designVariable = 'Density';
filterType = 'P1';

nsteps = 50;
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
maxiter = 100;
monitoring_interval = 1;
