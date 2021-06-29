%filename='Cantilever_tetrahedra';
%filename='Cantileverbeam_Tetrahedra_Linear_Structured_Fine';
%filename = 'Cantilever3D';
%filename = 'Cantilever3DLarge';
filename = 'Cantilever3DLargeFineFine';
%filename='Cantilever_quad_fine';
 ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
%cost = {'compliance','perimeter'};
cost = {'compliance'};
%weights = [1, 1];
weights = 1;
constraint = {'volumeConstraint'};
optimizerUnconstrained = 'SLERP'; 
incrementFactor = 1;
designVariable = 'LevelSet';

%optimizerUnconstrained = 'PROJECTED GRADIENT'; 
%designVariable = 'Density';

filterType = 'PDE';
optimizer = 'DualNestedInPrimal';
%optimizer = 'AlternatingPrimalDual';
line_search_initiator = 'INCREASING LAST STEP';
incremenFactor = 1.2;

nsteps = 85;
Vfrac_final = 0.15;
Perimeter_target = 1;
optimality_final = 1e-4;
constr_final = 1e-4;

Vfrac_initial      = 1;
optimality_initial = 1*1e-4;
constr_initial     = 1e-4;


TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

plotting = true;
printing = true;
monitoring = true;
monitoring_interval = 1;

maxiter = 8500;