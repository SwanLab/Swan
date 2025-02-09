filename='CantileverBeam_Quadrilateral_Bilinear_Structured';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'holes';
cost = {'compliance';'perimeter'};
weights = [1 0.1];
constraint = {'volumeConstraint'};
optimizer = 'SLERP'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';
constraints_case = 'INEQUALITY';
line_search_initiator = 'INCREASING LAST STEP';

nsteps = 10;
Vfrac_final = 0.2;
Perimeter_target=1;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.5;
optimality_initial = 1e-2;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;
