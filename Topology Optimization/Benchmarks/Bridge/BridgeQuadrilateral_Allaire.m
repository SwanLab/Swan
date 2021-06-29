filename='Bridge_Quadrilateral_Bilinear_Structured';
ptype = 'MACRO';
method = 'SIMP_P3'; % !! Instead of proportional to material density !!
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance','perimeter'};
weights = [1 0.02];
constraint = {'volumeConstraint'};
optimizer = 'HAMILTON-JACOBI'; incrementFactor = 1;
filterType = 'P1';
constraint_case = 'INEQUALITY';

HJiter0 = 1;
e2 = 50;
N_holes = [5 6];
R_holes = 0.7;
phase_holes = [pi pi];

nsteps = 10;
Vfrac_final = 0.1;
Perimeter_target=3.5;
optimality_final =5e-4;
constr_final =5e-4;

Vfrac_initial = 0.6;
optimality_initial = 5e-2;
constr_initial = 5e-2;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 0.3;
TOL.nu_minus = 0.3;

plotting = true;
showBC = false;
printing = false;
monitoring = true;
monitoring_interval = 1;