filename='test2d_micro';
%filename='RVE_Square_Triangle_FineFine';
%filename='RVE_Square_Triangle_FineFine';
%filename = 'MicroQuad';
ptype = 'MICRO';
method = 'SIMPALL';
%method = 'SIMP_P3';
materialType = 'ISOTROPIC';
%initial_case = 'circleInclusion';
%fracRadius = 0.95;
%fracRadius = 0.8;
%fracRadius = 0.4;
%initial_case = 'rectangleInclusion';
%widthH = 1.1;0.5;
%widthV = 0.5;
initial_case = 'holes';
N_holes = [8 8 4];
R_holes = 0.8;
phase_holes = [pi/2 pi/2 0];

%cost={'enforceCh_CCstar_L2'};
%cost={'chomog_alphabeta'};
%cost={'chomog_fraction'};
%weights = 1;
cost={'enforceCh_CCstar_L2','perimeterInterior'};
%weights=[1,1e-15];
weights=[1,1e-3];
constraint = {'volumeConstraint'};
constraint_case = 'INEQUALITY';
%constraint_case = 'EQUALITY';
%incrementFactor = 1;
designVariable = 'Density';
%designVariable = 'LevelSet';
%filterType = 'P1';
filterType = 'PDE';
%optimizer = 'IPOPT';
%optimizer = 'MMA';
optimizer = 'AlternatingPrimalDual';

%optimizer = 'DualNestedInPrimal';
optimizerUnconstrained = 'PROJECTED GRADIENT';
%optimizerUnconstrained = 'SLERP';

line_search_initiator = 'INCREASING LAST STEP';
incrementFactor = 2.95;
%incrementFactor = 1.5;


%selectiveC_Cstar = 'Composite';
%selectiveC_Cstar = 'Vfrac04b';
%selectiveC_Cstar = 'HorizontalRectangleInclusion';
%selectiveC_Cstar = 'NegativePoiss06';
%selectiveC_Cstar = 'nu_0_6';
%selectiveC_Cstar = 'Seba';
selectiveC_Cstar = 'Nu0_2';

nsteps = 6;
Vfrac_final = 0.99;
Perimeter_target=1;
optimality_final = 0.2*1e-3;
constr_final =1e-12;
%constr_final =1e-1;

Vfrac_initial = 0.99;
optimality_initial = 0.2*1e-3;
constr_initial = 1e-12;
%constr_initial = 1e-1;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

%Micro
epsilon_isotropy_initial=1e-1;
epsilon_isotropy_final = 1e-3;
micro.alpha = [1 0 0]';
micro.beta  = [0 -1 0]';

%epsilonPerFinal   = 2*m.computeMeanCellSize();
%epsilonPerInitial = 124*m.computeMeanCellSize();  

%nsteps = 10;

% For all tests
plotting = true;
printing = true;
monitoring = true;
monitoring_interval = 1;
maxiter = 500;%2000;