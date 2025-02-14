filename='Cantilever_triangle_fine';
 ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance'};
weights = 1;
constraint = {'volume','perimeter'};
optimizer = 'HAMILTON-JACOBI'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';
constraint_case = 'INEQUALITY';
line_search_initiator = 'STANDARD';

HJiter0 = 1;
e2 = 200;
N_holes = [1 2];
R_holes = 0.9;
phase_holes = [pi/2 pi/2];


nsteps = 15;
Vfrac_final = 0.4;
Perimeter_target = 1.2;
optimality_final = 1e-3;
constr_final =1e-3;

Vfrac_initial = 0.6;
optimality_initial = 1e-3;
constr_initial = 1e-1;
TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

monitoring_interval = 1;