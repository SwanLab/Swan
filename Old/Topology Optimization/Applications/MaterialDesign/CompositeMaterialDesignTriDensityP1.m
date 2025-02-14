filename='test2d_micro';
%filename='RVE_Square_Triangle_FineFine';
%filename = 'MicroQuad';
ptype = 'MICRO';
method = 'SIMPALL';
%method = 'SIMP_P3';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
cost={'chomog_alphabeta'};
weights=[1];
constraint = {'volumeConstraint'};
constraint_case = 'EQUALITY';
%incrementFactor = 1;
designVariable = 'Density';
%designVariable = 'LevelSet';
filterType = 'P1';
fracRadius = 0.4;
%optimizer = 'IPOPT';
optimizer = 'DualNestedInPrimal';
%optimizer = 'AlternatingPrimalDual';

optimizerUnconstrained = 'PROJECTED GRADIENT';
line_search_initiator = 'INCREASING LAST STEP';
incrementFactor = 2.95;

%optimizerUnconstrained = 'SLERP';


nsteps = 1;
Vfrac_final = 0.25;
Perimeter_target=1;
optimality_final = 0.2*1e-3;
constr_final =1e-12;

Vfrac_initial = 0.8;
optimality_initial = 0.2*1e-3;
constr_initial = 1e-12;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

%Micro
epsilon_isotropy_initial=1e-1;
epsilon_isotropy_final = 1e-3;
micro.alpha =[1 1 0]';
micro.beta =[1 1 0]';

%nsteps = 10;

% For all tests
plotting = true;
printing = true;
monitoring = true;
monitoring_interval = 1;
maxiter = 2000;