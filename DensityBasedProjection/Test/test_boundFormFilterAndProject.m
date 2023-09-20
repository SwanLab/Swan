filename = 'MeshReplicated';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'LinearBoundFunction'};
widthSquare = 0.01;
weights = [1];
constraint = {'ComplianceConstraintThreeFieldRhoE','ComplianceConstraintThreeFieldRhoI',...
    'ComplianceConstraintThreeFieldRhoD','VolumeConstraintRhoD'};
constraint_case = {'INEQUALITY','INEQUALITY','INEQUALITY','INEQUALITY'};
optimizerUnconstrained = 'PROJECTED GRADIENT';
optimizer = 'MMA';
incrementFactor = 1.2;
designVariable = 'Density&Bound';
filterType = 'Filter&Project';
nsteps = 1;
Vfrac_final = 0.5;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

% For all tests
plotting = false;
printing = false;
printing_physics = false;
monitoring = false;
monitoring_interval = 1;
maxiter = 5;