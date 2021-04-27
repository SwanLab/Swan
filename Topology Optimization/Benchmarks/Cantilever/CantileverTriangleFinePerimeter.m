filename='Cantilever_triangle_fine';
%filename='Cantilever_quad_fine';
 ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance','perimeterInterior'};
%cost = {'compliance','perimeter'};
%cost = {'compliance'};
weights = [1, 0.2];
%weights = 1;
constraint = {'volumeConstraint'};
optimizerUnconstrained = 'SLERP'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';
optimizer = 'DualNestedInPrimal';
%optimizer = 'AlternatingPrimalDual';


nsteps = 5;
Vfrac_final = 0.4;
Perimeter_target = 1;
optimality_final = 1e-4;
constr_final =1e-4;

Vfrac_initial = 0.4;
optimality_initial = 1e-4;
constr_initial = 1e-4;
TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

plotting = true;
printing = true;
monitoring = true;
monitoring_interval = 15;

maxiter = 300;