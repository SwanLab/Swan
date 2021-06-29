filename='CantileverSquare';
ptype = 'MACRO';
initial_case = 'full';
%cost = {'compliance'};
cost = {'stressNorm'};
weights = 1;
constraint = {'volumeConstraint'};
filterType = 'P1';

optimizer = 'DualNestedInPrimal';
optimizerUnconstrained = 'PROJECTED GRADIENT';
designVariable = 'MicroParams';
homegenizedVariablesComputer = 'ByVademecum';
vademecumFileName = 'Rectangle';%'SmoothRectangle';
%vademecumFileName = 'SmoothRectangle';

nsteps = 1;
maxiter = 500;

Vfrac_initial = 0.3;
Vfrac_final = 0.3;

optimality_initial = 1e-5;
optimality_final = 1e-5;
constr_initial = 1e-3;
constr_final = 1e-3;

printing = false;
monitoring_interval = 1;

ub = 0.989;
lb = 0.011;
rate = 0.5;




