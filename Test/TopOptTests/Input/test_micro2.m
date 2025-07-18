filename = 'test2d_micro';
a.fileName = filename;
gid = FemDataContainer(a);
ptype = 'MICRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
geomFunSettings.type = 'CircleInclusion';
geomFunSettings.xCoorCenter = 0.5;
geomFunSettings.yCoorCenter = 0.5;
geomFunSettings.radius      = 0.2;
cost = {'chomog_alphabeta'};
weights = [1];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
target = 0.5;
optimizerUnconstrained = 'SLERP'; 
optimizer = 'AlternatingPrimalDual';
designVariable = 'LevelSet';
filterCostType = {'P1'};
filterConstraintType = {[]};
filterCostSettings = {[]};
filterConstraintSettings = {[]};
micro.alpha =[1,0; 0,1];
micro.beta =[1,0; 0,1];
plotting = false;
printing = false;
monitoring = false;
maxiter = 5;