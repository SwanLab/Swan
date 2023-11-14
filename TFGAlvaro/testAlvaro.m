% ONLY CHANGE PARAMETERS WITH COMMENTS

filename = 'CantileverAlvaro'; % Diff meshes
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
optimizerUnconstrained = 'PROJECTED GRADIENT'; % PROJECTED GRADIENT (Density)/ SLERP (LevelSet)
optimizer = 'NullSpace';
incrementFactor = 1.2;
designVariable = 'Density'; % Density / LevelSet
filterType = 'P1';
nsteps = 1;
Vfrac_final = 0.3; % Final volume fraction
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

plotting = false;
printing = false; % Maybe true to obtain evolution(?)
printing_physics = false;
monitoring = true;
monitoring_interval = 1;
maxiter = 300; % Maximum number of iterations