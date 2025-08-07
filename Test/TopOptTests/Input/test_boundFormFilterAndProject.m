filename = 'anisoCantilever';
a.fileName = filename;
gid = FemDataContainer(a);
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
geomFunSettings.type = 'Full';
cost = {'LinearBoundFunction'};
weights = [1];
constraint = {'ComplianceConstraintBound','ComplianceConstraintBound',...
    'ComplianceConstraintBound','VolumeConstraintBound'};
constraint_case = {'INEQUALITY','INEQUALITY','INEQUALITY','INEQUALITY'};
target = [NaN,NaN,NaN,0.5];
optimizerUnconstrained = 'PROJECTED GRADIENT';
optimizer = 'MMA';
designVariable = 'DensityAndBound';
filterCostType = {[]};
filterConstraintType = {'FilterAndProject','FilterAdjointAndProject',...
    'FilterAndProject','FilterAdjointAndProject','FilterAndProject',...
    'FilterAdjointAndProject',[]};
filterCostSettings = {[]};
f1.filterStep = 'PDE';
f1.eta = 0.75;
f1.beta = 1;
f2.filterStep = 'PDE';
f2.eta = 0.5;
f2.beta = 1;
f3.filterStep = 'PDE';
f3.eta = 0.25;
f3.beta = 1;
filterConstraintSettings = {f1,f1,f2,f2,f3,f3,[]};
plotting = false;
printing = false;
monitoring = false;
maxiter = 20;