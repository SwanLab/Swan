filename = 'Cantileverbeam_Quadrilateral_Bilinear';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initialCase = 'Full';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
target = 0.1;
optimizerUnconstrained = 'SLERP';
optimizer = 'NullSpace';
designVariable = 'LevelSet';
filterType = 'P1';

E1  = 1;
E0  = 1e-3;
nu1 = 1/3;
nu0 = 1/3;

plotting = false;
printing = false;
monitoring = false;
maxiter = 3;