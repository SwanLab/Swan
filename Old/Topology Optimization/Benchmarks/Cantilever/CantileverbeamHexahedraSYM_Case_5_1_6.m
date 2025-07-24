filename='Cantileverbeam_Hexahedra_Linear_Structured_SYM';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance'};
weights = [1];
constraint = {'volume'};
optimizer = 'HAMILTON-JACOBI'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';
constraint_case = 'INEQUALITY';

nsteps = 10;
Vfrac_final = 0.1;
Perimeter_target=3.5;
optimality_final =1e-3;
constr_final =1e-3;

BCscale_factor = 0.3;
rotation_per_it = 20;
HJiter0 = 1;
e2 = 100;
N_holes = [6 3 3];
R_holes = 0.4;
phase_holes = [0 0 0];

Vfrac_initial = 0.7;
optimality_initial = 5e-2;
constr_initial = 5e-2;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;

% maxiter = 1;
TOL.nu_plus = 0.3;
TOL.nu_minus = 0.3;