filename = 'Cantileverbeam_Quadrilateral_Bilinear';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
geomFunSettings.type = 'RectangleInclusion';
geomFunSettings.xSide = 0.8;
geomFunSettings.ySide = 0.4;
geomFunSettings.xCoorCenter = 1;
geomFunSettings.yCoorCenter = 0.5;
cost = {'compliance';'perimeter'};
weights = [1 0.1];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
target = 0.3;
optimizer = 'SLERP';
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';

shFuncParamsName = 'paramsTestCantilever3';

nsteps = 1;
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
plotting = false;
printing = false;
printing_physics = false;
monitoring = false;
maxiter = 3;
