filename='RVE_Square_Triangle';
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
cost={'enforceCh_CCstar_L2','perimeter'};%enforceCh_CCstar_L2
weights=[1, 0.1];
constraint = {'volumeConstraint'};
optimizer = 'MMA'; incrementFactor = 1;%mma CON 1E-4
filterType = 'P1';

nsteps = 5;
Vfrac_final = 0.6;
Perimeter_target=1;
optimality_final =1e-4;
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
selectiveC_Cstar = 'Vfrac06';