filename = 'Gripping_triangle_coarse';
a.fileName = filename;
gid = FemDataContainer(a);
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
geomFunSettings.type = 'Full';
cost = {'nonadjointCompliance'};
weights = [1];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
target = 0.5;
optimizerUnconstrained = 'SLERP';
optimizer = 'NullSpace';
designVariable = 'LevelSet';
filterCostType = {'P1'};
filterConstraintType = {[]};
filterCostSettings = {[]};
filterConstraintSettings = {[]};
plotting = false;
printing = false;
monitoring = false;
maxiter = 5;