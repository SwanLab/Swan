filename = 'CantileverSquareSmall';
ptype = 'MACRO';
initial_case = 'full';
cost = {'stressNorm'};
weights = 1;
constraint = {'volumeConstraint'};
filterType = 'P1';
constraint_case = 'EQUALITY';

m1 = 0.3;
m2 = 0.3;

stressNormExponent_initial = 2;
stressNormExponent_final = 16;

Vfrac_initial = 0.3;
optimality_initial = 1e-3;
constr_initial = 1e-3;

Vfrac_final = 0.3;
optimality_final = 1e-5;
constr_final = 1e-3;

optimizer = 'DualNestedInPrimal';
optimizerUnconstrained = 'PROJECTED GRADIENT';
designVariable = 'MicroParams';
ub = 0.989;
lb = 0.011;
kfrac = 2;
nsteps = 8;
homegenizedVariablesComputer = 'ByVademecum';
vademecumFileName = 'SuperEllipseQOptAnalytic';


plotting = false;
printing = false;
printing_physics = false;
monitoring = false;
monitoring_interval = 1;
maxiter = 16;
