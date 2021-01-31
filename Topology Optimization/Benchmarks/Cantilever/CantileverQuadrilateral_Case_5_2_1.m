filename='CantileverBeam_Quadrilateral_Bilinear_Structured';
ptype = 'MACRO';
method = 'SIMPALL'; % !! Instead of proportional to material density !!
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance','perimeter'};
weights = [1 0.1];
constraint = {'volumeConstraint'};
optimizer = 'HAMILTON-JACOBI'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';
constraint_case = 'INEQUALITY';
line_search_initiator = 'INCREASING LAST STEP';

nsteps = 10;
Vfrac_final = 0.1;
Perimeter_target=3.5;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.4;
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