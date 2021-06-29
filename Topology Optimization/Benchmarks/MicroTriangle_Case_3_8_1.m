filename='RVE_Square_Triangle';
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
cost={'volume','perimeter'};
weights=[1 0.1];
constraint = {'enforceCh_CCstar_eq'};%inf & equality
optimizer = 'MMA'; incrementFactor = 1;
filterType = 'P1';

nsteps = 1;
Vfrac_final = 1;
Perimeter_target=1;
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
selectiveC_Cstar = [1,1,1;
    1,1,1;
    1,1,1]; % 1 to select the component



