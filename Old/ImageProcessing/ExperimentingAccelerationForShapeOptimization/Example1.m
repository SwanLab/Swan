filename='CantileverBeam_Triangle_Linear_Fine';
ptype = 'MACRO';
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'full';
cost = {'compliance'};
weights = 1;
%cost = {'compliance'};
%weights = [1];
constraint = {'volumeConstraint'};
%optimizer = 'DualNestedInPrimal';
optimizer = 'AlternatingPrimalDual';

optimizerUnconstrained = 'PROJECTED GRADIENT'; 
incrementFactor = 1;
designVariable = 'Density';
filterType = 'P1';

nsteps = 1;
Vfrac_final = 0.4;
optimality_final =1e-3;
constr_final = 1e-3;

Vfrac_initial = 1;
optimality_initial = 1e-2;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

plotting = false;
printing = false;
monitoring = false;

maxiter = 200;

monitoring_interval = 1;