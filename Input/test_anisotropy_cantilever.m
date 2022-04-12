% filename = 'Cantileverbeam_Quadrilateral_Bilinear';
filename = 'CantileverBeam_Triangle_Linear';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance';'perimeter'};
weights = [1 0.1];
constraint = {'volumeConstraint'};
optimizer = 'SLERP'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'AnisotropicPDE';

shFuncParamsName = 'paramsTestCantilever3';

nsteps = 5;
Vfrac_final = 0.3;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-3;
constr_initial = 1e-3;
Perimeter_target = 1;

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
maxiter = 150;
