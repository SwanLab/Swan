filename='Bridge_Quadrilateral_Bilinear_Structured';
ptype = 'MACRO';
method = 'SIMP_P3'; % !! Instead of proportional to material density !!
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance'};
weights = [1];
constraint = {'volume'};
optimizer = 'HAMILTON-JACOBI'; kappaMultiplier = 1;
filterType = 'PDE';

nsteps = 1;
Vfrac_final = 0.5;
Perimeter_target=3.5;
optimality_final =1e-3;
constr_final =1e-3;
HJiter0 = 30;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 0.3;
TOL.nu_minus = 0.3;

plotting = 1;
printing = 0;
monitoring = 1;
monitoring_interval = 1;