filename = 'Cantileverbeam_Tetrahedra_Linear_Structured_Fine';%'Cantileverbeam_Quadrilateral_Bilinear';%CantileverArnau2 %Cantilever';%'CantileverBeam_Triangle_Linear';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';%'holes';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
optimizerUnconstrained = 'PROJECTED GRADIENT';%'SLERP';%'PROJECTED GRADIENT';%'PROJECTED GRADIENT'; 
optimizer = 'NullSpace';%'DualNestedInPrimal';'DualNestedInPrimal';%'AlternatingPrimalDual';%'AlternatingPrimalDual';
incrementFactor = 1;
designVariable = 'Density';%'Density';'LevelSet'
filterType = 'P1';
nsteps = 1;
Vfrac_final = 0.4;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;
Perimeter_target = 5;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

% For all tests
plotting = true;
printing = false;
printing_physics = false;
monitoring = true;
monitoring_interval = 5;
maxiter = 300;