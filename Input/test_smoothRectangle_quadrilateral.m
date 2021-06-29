filename='Circle_Quadrilateral_Linear_Structured_32';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circle';
fracRadius = 1.2;
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
optimizer = 'HAMILTON-JACOBI';
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';

nsteps = 1;
Vfrac_final = 1;
Perimeter_target=3.5;
optimality_final =1e-5;
constr_final =1e-5;

BCscale_factor = 0.3;
HJiter0 = 1;
e2 = 5;

Vfrac_initial = 1;
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

plotting = 0;
printing = 0;
monitoring = 0;
monitoring_interval = 0;