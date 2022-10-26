filename = 'Cantilever';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
constraint_case = 'EQUALITY';
optimizerUnconstrained = 'SLERP'; 
optimizer = 'AlternatingPrimalDual';
incrementFactor = 1.05; 
designVariable = 'LevelSet';
filterType = 'PDE';
line_search_initiator = 'INCREASING LAST STEP';


nsteps = 1;
Vfrac_final = 0.3;
optimality_final =1e-3;
constr_final = 1e-5;

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

plotting = false;
printing = false;
printing_physics = false;
monitoring = false;
monitoring_interval = 10;
maxiter = 15;
