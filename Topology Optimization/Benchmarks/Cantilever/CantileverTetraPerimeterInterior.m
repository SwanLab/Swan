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
cost = {'compliance','perimeterInterior'};
%cost = {'compliance'};
weights = [1, 1];
%weights = 1;
constraint = {'volumeConstraint'};
optimizerUnconstrained = 'SLERP'; 
incrementFactor = 1;
designVariable = 'LevelSet';

%optimizerUnconstrained = 'PROJECTED GRADIENT'; 
%designVariable = 'Density';

filterType = 'PDE';
optimizer = 'DualNestedInPrimal';
%optimizer = 'AlternatingPrimalDual';


nsteps = 16;
Vfrac_final = 0.15;
Perimeter_target = 1;
optimality_final = 1e-3;
constr_final = 1e-2;

Vfrac_initial      = 0.15;
optimality_initial = 0.5*1e-2;
constr_initial     = 1e-2;


TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

plotting = true;
printing = false;
monitoring = true;
monitoring_interval = 1;

maxiter = 160;