%filename = 'CantileverBeam_Triangle_Linear_FineFine';
%filename = 'CantileverBeam_Triangle_Linear_Fine';
%filename = 'Cantilever_quad_fine';
%filename = 'CantileverLargeCoarse';
filename = 'CantileverLargeFine';
ptype = 'MACRO';
initial_case = 'full';
cost = {'compliance','perimeterInterior'};
weights = [1 0.1];
%cost = {'compliance'};
%cost = {'stressNorm'};
%weights = [1];
constraint = {'volumeConstraint'};
optimizer = 'DualNestedInPrimal';
%optimizer = 'AlternatingPrimalDual';

optimizerUnconstrained = 'SLERP'; 
designVariable = 'LevelSet';
method = 'SIMPALL';
materialType = 'ISOTROPIC';

% optimizerUnconstrained = 'PROJECTED GRADIENT';
% designVariable = 'Density';
% method = 'SIMPALL';
% materialType = 'ISOTROPIC';

% optimizerUnconstrained = 'PROJECTED GRADIENT';
% designVariable = 'MicroParams';
% ub = 0.989;
% lb = 0.011;
% homegenizedVariablesComputer = 'ByVademecum';
% vademecumFileName = 'SuperEllipseQOptAnalytic';
% m1 = 0.0101;
% m2 = 0.0101;



%filterType = 'P1';
filterType = 'PDE';

nsteps = 20;
Vfrac_final = 0.4;
optimality_final =1e-4;
constr_final = 1e-4;

Vfrac_initial = 1;
optimality_initial = 1e-2;
constr_initial = 1e-3;

TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;

plotting = true;
printing = true;
monitoring = true;

maxiter = 200;

monitoring_interval = 1;