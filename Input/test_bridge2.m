filename='Bridge_Tetrahedra_Coarse';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance'};
weights = 1;
constraint = {'volumeConstraint'};
optimizer = 'SLERP'; 
incrementFactor = 1;
filterType = 'P1';
designVariable = 'LevelSet';
constraint_case = 'INEQUALITY';

nsteps = 1;
Vfrac_final = 0.1;
Perimeter_target=3.5;
optimality_final =1e-5;
constr_final =1e-5;

BCscale_factor = 0.3;
HJiter0 = 1;
e2 = 1;
N_holes = [6 3 3];
R_holes = 0.9;
phase_holes = [0 0 0];

Vfrac_initial = 0.3;
optimality_initial = 5e-2;
constr_initial = 5e-2;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;

maxiter = 3;
TOL.nu_plus = 0.3;
TOL.nu_minus = 0.3;

plotting = false;
printing = false;
monitoring = false;
monitoring_interval = 1;