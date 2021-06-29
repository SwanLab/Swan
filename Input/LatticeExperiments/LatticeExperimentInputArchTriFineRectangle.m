filename = 'ArchTriFine';%'CantileverSquareNewFine';%'LshapeTriFine';%'CantileverSquareNew';'ArchTriFine';%'CantileverSquare';%'ArchTriFine';%'CantileverSquare';%'Lshape';%'LshapeTriFine';%'CantileverSquareSmall';%'';%'Lshape';'CantileverSquare';'Lshape';%'CantileverSquare';%'LshapeFine';%'Bridge_quad_coarse';%'Arch_quad_coarse';%'BridgeCool_Quadrilateral_Bilinear_Structured_Coarse';%'Bridge';%'CantileverSquareSmall';
ptype = 'MACRO';
initial_case = 'given';
m1 = 0.0101;
m2 = 0.0101;
%cost = {'compliance'};
cost = {'stressNorm'};
%cost = {'stressNorm','compliance'};
%weights = [0.55,0.45];
weights = 1;
constraint = {'volumeConstraint'};
filterType = 'P1';
constraint_case = 'EQUALITY';

Vfrac_initial = 0.3;
optimality_initial = 1e-5;
constr_initial = 1e-5;

Vfrac_final = 0.3;
optimality_final = 1e-5;
constr_final = 1e-5;

stressNormExponent_initial = 2;
stressNormExponent_final = 16;

optimizer = 'DualNestedInPrimal';
optimizerUnconstrained = 'PROJECTED GRADIENT';

designVariable = 'MicroParams';
ub = 0.989;
lb = 0.011;
homegenizedVariablesComputer = 'ByVademecum';
% 
vademecumFileName = 'SuperEllipseQMax';
%vademecumFileName = 'SuperEllipseQ2';
%vademecumFileName = 'SuperEllipseQOptAnalytic';
% 
% designVariable = 'Density';
% homegenizedVariablesComputer = 'ByInterpolation';
% method = 'SIMPALL';
% materialType = 'ISOTROPIC';

kfrac = 2;
nsteps = 10;

plotting = false;
printing = true;
monitoring = true;
monitoring_interval = 1;
maxiter = 500;


