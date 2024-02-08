filename = 'CantileverBeam_Triangle_Linear';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initialCase = 'Full';
cost = {'compliance', 'perimeter'};
weights = [1 0.1];
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
target = 0.3;
optimizerUnconstrained = 'PROJECTED GRADIENT'; 
optimizer = 'AlternatingPrimalDual';
designVariable = 'Density';
filterCostType = {'P1','PDE'};
filterConstraintType = {[]};

E1  = 1;
E0  = 1e-3;
nu1 = 1/3;
nu0 = 1/3;

plotting = false;
printing = false;
monitoring = false;
maxiter = 15;