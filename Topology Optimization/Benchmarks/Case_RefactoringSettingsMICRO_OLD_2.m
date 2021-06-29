filename = 'test2d_micro';
ptype = 'MICRO';
method = 'SIMP_P3';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
cost = {'chomog_fraction';'perimeter'};
weights = [1 0.1];
constraint = {'volumeConstraint'};
constraint_case = 'INEQUALITY';
optimizer = 'MMA'; 
incrementFactor = 1;
designVariable = 'Density';
filterType = 'P1';

nsteps = 1;
Vfrac_final = 0.5;
optimality_final =1e-3;
constr_final =1e-3;

Vfrac_initial = 0.5;
optimality_initial = 1e-2;
constr_initial = 1e-2;
% Perimeter_target = 5;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

 %Micro
target_parameters.epsilon_isotropy = 1e-3;%%%%%%
epsilon_isotropy_final=target_parameters.epsilon_isotropy;
epsilon_isotropy_initial=1e-1;
micro.alpha =[1 0 0]';%[1 0 0]'
micro.beta =[0 -1 0]';%[0 -1 0]'

% For all tests
plotting = true;
printing = false;
printing_physics = false;
monitoring = true;

isOld = false;
