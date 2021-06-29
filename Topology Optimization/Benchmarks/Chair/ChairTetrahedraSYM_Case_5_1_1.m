filename='Chair_Tetrahedra_SYM';
 ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance'};
weights = [1];
constraint = {'volume'};
optimizer = 'HAMILTON-JACOBI'; incrementFactor = 1;
filterType = 'P1';
line_search_initiator = 'STANDARD';

HJiter0 = 1;
e2 = 100;
N_holes = [24 11 11];
R_holes = 0.4;
phase_holes = [0 0 0];
nsteps = 5;

Vfrac_final = 0.2;
Perimeter_target = 1;
optimality_final = 1e-3;
constr_final =1e-3;

Vfrac_initial = 0.4;
optimality_initial = 5e-2;
constr_initial = 5e-2;
TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;