filename='Bridge_Quadrilateral_Bilinear_Structured';
ptype = 'MACRO';
method = 'SIMPALL'; % !! Instead of proportional to material density !!
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint','perimeterConstraint'};
optimizer = 'HAMILTON-JACOBI'; incrementFactor = 1;
filterType = 'P1';

nsteps = 1;
Vfrac_final = 0.5;
Perimeter_target=3.5;
optimality_final =1e-3;
constr_final =1e-3;
HJiter0 = 30;
N_holes = [5 6];
R_holes = 0.7;
phase_holes = [0 pi/2];

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