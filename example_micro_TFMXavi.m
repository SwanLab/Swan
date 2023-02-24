% INPUTS:
%   General:
filename = 'instantSquare'; % MESH
ptype = 'MICRO'; % MICRO (unit cell) or MACRO (structure)
cost = {'chomog_alphabeta'};
weights = [1];
designVariable = 'Density'; % Density or LevelSet
fracRadius = 0.5; % Only needed for MICRO
Vfrac_final = 0.5;
%   Micro:
micro.alpha =[0 0 1]';
micro.beta =[0 0 1]';
maxiter = 1000;


% TopOpt stuff (may not change):
method = 'SIMPALL';
materialType = 'ISOTROPIC';
initial_case = 'circleInclusion';
constraint = {'volumeConstraint'};
constraint_case = {'EQUALITY'};
optimizer = 'DualNestedInPrimal';
switch designVariable
    case 'Density'
        optimizerUnconstrained = 'PROJECTED GRADIENT';
    case 'LevelSet'
        optimizerUnconstrained = 'SLERP';
end
incrementFactor = 1;
filterType = 'P1';
nsteps = 1;
optimality_final =1e-3;
constr_final =1e-3;
Vfrac_initial = 1-pi*fracRadius^2;
optimality_initial = 1e-3;
constr_initial = 1e-3;
Perimeter_target = 5;
TOL.rho_plus = 1;
TOL.rho_minus = 0;
TOL.E_plus = 1;
TOL.E_minus = 1e-3;
TOL.nu_plus = 1/3;
TOL.nu_minus = 1/3;
target_parameters.epsilon_isotropy = 1e-3;
epsilon_isotropy_final=target_parameters.epsilon_isotropy;
epsilon_isotropy_initial=1e-1;
plotting = true;
printing = false;
printing_physics = false;
monitoring = true;