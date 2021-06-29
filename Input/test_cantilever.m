filename = 'CantileverBeam_Triangle_Linear';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance','perimeter'};
weights = [1 0.1];
constraint = {'volumeConstraint'};
constraint_case = 'EQUALITY';
optimizerUnconstrained = 'PROJECTED GRADIENT'; 
optimizer = 'AlternatingPrimalDual';
incrementFactor = 1.05; 
designVariable = 'Density';
filterType = 'P1';
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
monitoring_interval = 1;
maxiter = 15;
