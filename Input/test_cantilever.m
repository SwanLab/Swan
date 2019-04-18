filename = 'Tests_Triangle_Linear';
ptype = 'MACRO';
method = 'SIMP_Adaptative';
materialType = 'ISOTROPIC';
initial_case = 'square';
widthSquare = 0.4;
cost = {'compliance';'perimeter'};
weights = [1 0.1];
constraint = {'volumeConstraint'};
constraint_case = 'INEQUALITY';
optimizer = 'PROJECTED GRADIENT'; 
kappaMultiplier = 1; 
designVariable = 'Density';
filterType = 'PDE';

shFuncParamsName = 'paramsTestCantilever';

nsteps = 1;
Vfrac_final = 1;
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

plotting = false;
printing = false;
printing_physics = false;
monitoring = false;
maxiter = 3;
