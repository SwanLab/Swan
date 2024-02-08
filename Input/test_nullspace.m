filename = 'Cantileverbeam_Quadrilateral_Bilinear';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
geomFunSettings.type = 'Full';
cost = {'compliance'};
weights = [1];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
target = 0.1;
optimizerUnconstrained = 'SLERP';
optimizer = 'NullSpace';
designVariable = 'LevelSet';
filterCostType = {'P1'};
filterConstraintType = {[]};

E1  = 1;
E0  = 1e-3;
nu1 = 1/3;
nu0 = 1/3;

plotting = false;
printing = false;
monitoring = false;
maxiter = 3;