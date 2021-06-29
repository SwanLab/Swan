filename='RVE_Square_Triangle';
ptype = 'MICRO';
method = 'SIMP_P3';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
cost={'chomog_fraction';'perimeter'};
weights=[1 0.1];
constraint = {'volumeConstraint'};
optimizer = 'MMA'; incrementFactor = 1;
filterType = 'P1';

nsteps = 1;
Vfrac_final = 0.5;
Perimeter_target=3.5;
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

%Micro
epsilon_isotropy_initial=1e-1;
epsilon_isotropy_final = 1e-3;
micro.alpha =[1 0 0]';%[1 0 0]'
micro.beta =[0 -1 0]';%[0 -1 0]'

% printing = true;
% maxiter = 15;


