%filename = 'CantileverSquareSYmmetricMesh';
%filename = 'CantileverSquareMedium';
%'CantileverSquareSYmmetricMesh';
%filename = 'CantileverSquare';
filename = 'CantileverSquareSmall';
%filename = 'CantileverSquareNew';
%'CantileverSquareNewFine';
%'LshapeTriFine';
%'Lshape';
%'LshapeFine';
%'ArchTriFine';
%'Arch_quad_coarse';
%'Bridge_quad_coarse';
%'BridgeCool_Quadrilateral_Bilinear_Structured_Coarse';
%'Bridge';


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
filterType = 'PDE';
constraint_case = 'EQUALITY';

Vfrac_initial = 0.3;
optimality_initial = 1e-5;
constr_initial = 1e-5;

Vfrac_final = 0.3;
optimality_final = 1e-5;
constr_final = 1e-5;

stressNormExponent_initial = 2;
stressNormExponent_final = 20;

optimizer = 'DualNestedInPrimal';
%optimizer = 'AlternatingPrimalDual';
optimizerUnconstrained = 'PROJECTED GRADIENT';

designVariable = 'MicroParams';
ub = 0.989;
lb = 0.011;
homegenizedVariablesComputer = 'ByVademecum';
% % 
vademecumFileName = 'SuperEllipseQMax';
%vademecumFileName = 'SuperEllipseQ2';
%vademecumFileName = 'SuperEllipseQOptAnalytic';
% 


%designVariable = 'LevelSet';
%homegenizedVariablesComputer = 'ByInterpolation';
%method = 'SIMPALL';
%materialType = 'ISOTROPIC';
%initial_case = 'full';
%optimizerUnconstrained = 'SLERP';


% designVariable = 'Density';
% homegenizedVariablesComputer = 'ByInterpolation';
% method = 'SIMPALL';
% materialType = 'ISOTROPIC';
% initial_case = 'full';



kfrac = 2;
nsteps = 10;

plotting = true;
printing = true;
monitoring = true;
monitoring_interval = 1;
maxiter = 5000;


