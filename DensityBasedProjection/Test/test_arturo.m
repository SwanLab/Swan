filename = 'Gripping25K';
%Gripping
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'LinearBoundFunction'};
widthSquare = 0.01;
weights = [1];
%cost = {'nonadjoint_compliance'};
constraint = {'ComplianceConstraintThreeFieldRhoE','ComplianceConstraintThreeFieldRhoI',...
    'ComplianceConstraintThreeFieldRhoD','VolumeConstraintRhoD'};
%constraint = {'volumeConstraint'};
%constraint_case = {'INEQUALITY'};
constraint_case = {'INEQUALITY','INEQUALITY','INEQUALITY','INEQUALITY'};
optimizerUnconstrained = 'PROJECTED GRADIENT';
optimizer = 'MMA';
incrementFactor = 1.2;
designVariable = 'Density&Bound';
%designVariable = 'Density';
filterType = 'Filter&Project';
%filterType = 'PDE';
nsteps = 1;
Vfrac_final = 0.35;
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
plotting = true;
printing = false;
printing_physics = false;
monitoring = true;
monitoring_interval = 1;
maxiter = 1500;