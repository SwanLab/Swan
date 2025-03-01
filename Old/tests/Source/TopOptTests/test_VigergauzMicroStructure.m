filename ='test2d_micro';%'RVE_Square_Triangle_FineFine';
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
cost={'chomog_alphabeta','perimeterConstraint'};
weights=[1 0.1];
constraint = {'volumeConstraint'};
optimizer = 'SLERP'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';
initial_case = 'Vigdergauz';
vigdergauzType = 'VolumeAndRatio';
superEllipseRatio = 0.9;
volumeMicro = 0.61;


nsteps = 1;
Vfrac_final = 0.5;
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
micro.alpha =[1 1 0]';
micro.beta =[1 1 0]';

% For all tests
plotting = false;
printing = false;
monitoring = false;
maxiter = 0;
