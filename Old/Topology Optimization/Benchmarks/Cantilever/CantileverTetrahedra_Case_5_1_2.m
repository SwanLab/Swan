filename='Cantilever_tetrahedra_coarse';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
optimizer = 'HAMILTON-JACOBI'; 
incrementFactor = 1;
designVariable = 'Density';
filterType = 'P1';
constraint_case = 'INEQUALITY';
line_search_initiator = 'INCREASING LAST STEP';

nsteps = 3;
Vfrac_final = 0.5;
Perimeter_target=3.5;
optimality_final =1e-5;
constr_final =1e-5;

BCscale_factor = 0.7;
HJiter0 = 1;
e2 = 1e1;
N_holes = [3 6 3];
R_holes = 0.3;
phase_holes = [0 0 0];

Vfrac_initial = 0.8;
optimality_initial = 5e-2;
constr_initial = 5e-2;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;

% maxiter = 10;
% maxiter = 1;
TOL.nu_plus = 0.3;
TOL.nu_minus = 0.3;

plotting = 1;
printing = 1;
monitoring = 1;
monitoring_interval = 1;