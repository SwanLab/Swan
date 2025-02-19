filename='Cantilever_hexahedra_coarse';
 ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance'};
weights = [1, 0.1];
constraint = {'volumeConstraint'};
optimizer = 'PROJECTED GRADIENT'; 
incrementFactor = 1;
designVariable = 'Density';
filterType = 'P1';

nsteps = 5;
Vfrac_final = 0.05;
Perimeter_target = 1;
optimality_final = 1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-2;
constr_initial = 1e-3;
TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;