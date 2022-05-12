filename = 'HexagonTest';
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion'; %% Parámetro a ir variando
cost = {'enforceCh_CCstar_L2';'perimeter'};
weights = [1 1e-3];
constraint = {'volumeConstraint'};
constraint_case = 'INEQUALITY';
optimizer = 'MMA'; 

incrementFactor = 1;
designVariable = 'Density';
filterType = 'P1';
fracRadius = 0.5;

nsteps = 1;
Vfrac_final = 0.5;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.5;
optimality_initial = 1e-3;
constr_initial = 1e-3;
Perimeter_target = 5;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

 %Micro
epsilon_isotropy_initial=1e-1;
epsilon_isotropy_final = 1e-3;
selectiveC_Cstar = 'IsotropyHexagon';

% For all tests
plotting = true;
printing = false;
printing_physics = false;
monitoring = true;
maxiter = 50;
