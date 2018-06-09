filename='Gripping_triangle_coarse';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'nonadjoint_compliance'};
weights = [1, 0.1];
constraint = {'volume'};
constraint_case = 'INEQUALITY';
optimizer = 'PROJECTED GRADIENT'; kappaMultiplier = 1;
filterType = 'P1';

nsteps = 1;
Vfrac_final = 1;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;
