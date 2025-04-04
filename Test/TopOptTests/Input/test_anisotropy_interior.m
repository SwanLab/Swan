filename = 'anisoCantilever';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
geomFunSettings.type = 'CircleInclusion';
geomFunSettings.xCoorCenter = 1;
geomFunSettings.yCoorCenter = 0.5;
geomFunSettings.radius      = 0.2;
cost = {'compliance','perimeter'};
weights = [1,0.025];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
target = 0.6;
optimizerUnconstrained = 'SLERP';
optimizer = 'NullSpace';
designVariable = 'LevelSet';
filterCostType = {'P1','PDE'};
filterConstraintType = {[]};
f1 = [];
scaleAngle = 60;
overhangAngle = 90;
tu = tand(scaleAngle);
f2.boundaryType = 'Neumann';
f2.metric       = 'Anisotropy';
f2.CAnisotropic = [tu,0;0,1/tu];
f2.aniAlphaDeg  = overhangAngle;
f3 = [];
filterCostSettings = {f1,f2};
filterConstraintSettings = {f3};
plotting = false;
printing = false;
monitoring = false;
maxiter = 5;