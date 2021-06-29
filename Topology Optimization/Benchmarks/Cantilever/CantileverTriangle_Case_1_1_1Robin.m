filename='CantileverBeam_Triangle_Linear_FineFine';
%filename = 'Cantilever_quad_fine';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance','perimeter'};
weights = [1 1];
constraint = {'volumeConstraint'};
optimizer = 'DualNestedInPrimal';
%weights = [1];
%constraint = {'volumeConstraint','perimeterConstraint'};
%optimizer = 'AlternatingPrimalDual';


optimizerUnconstrained = 'SLERP'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';

nsteps = 50;
Vfrac_final = 0.3;
optimality_final =1e-3;
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

epsilonPerInitial = 0.1;

maxiter = 1500;
