filename = 'inverter'; % inverter
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'nonadjoint_compliance'};
weights = [1];
constraint = {'volumeConstraint'};
incrementFactor = 1.2;
optimizerUnconstrained = 'SLERP';
optimizer = 'NullSpace';
designVariable = 'LevelSet';
filterType = 'PDE';
constraint_case = {'EQUALITY'};

nsteps = 1;
Vfrac_final = 0.18;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.18;
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
maxiter = 1000;