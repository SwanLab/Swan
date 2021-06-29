filename='SquareMacroStrucutred';
%filename='SquareDomainQuad';
ptype = 'MACRO';
initial_case = 'rectangleInclusion';
v = 0.4;
widthH = sqrt(1-v);
widthV = sqrt(1-v);
%cost = {'perimeter'};
cost = {'perimeterInterior'};
weights = 1;
constraint = {'volumeConstraint'};
optimizerUnconstrained = 'SLERP'; 
incrementFactor = 1;
designVariable = 'LevelSet';
filterType = 'P1';
optimizer = 'DualNestedInPrimal';

nsteps = 5;
Vfrac_initial = v;
Vfrac_final = v;
optimality_final = 1e-4;
constr_final =1e-6;

optimality_initial = 1e-4;
constr_initial = 1e-6;


plotting = true;
printing = false;
monitoring = true;
monitoring_interval = 1;

maxiter = 100;