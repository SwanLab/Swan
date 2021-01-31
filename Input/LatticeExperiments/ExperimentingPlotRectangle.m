%filename = 'CantileverSquareSYmmetricMesh';
%filename = 'CantileverSquareMedium';
%'CantileverSquareSYmmetricMesh';
%filename = 'CantileverSquare';
filename = 'CantileverSquareSmall';
%filename = 'CantileverSquareNew';
%'CantileverSquareNewFine';
%filename = 'Cantilever_quad_coarse';
%filename = 'LshapeTriFine';
%filename = 'LshapeTri';
%filename = 'LshapeTriSmall';
%filename = 'Lshape';
%filename = 'LshapeFine';
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
%filterType = 'P1';
constraint_case = 'EQUALITY';

Vfrac_initial = 0.3;
optimality_initial = 1e-3;
constr_initial = 1e-5;

Vfrac_final = 0.3;
optimality_final = 1e-3;
constr_final = 1e-5;

stressNormExponent_initial = 2;
stressNormExponent_final = 16;
% 
optimizer = 'DualNestedInPrimal';
%optimizer = 'AlternatingPrimalDual';
optimizerUnconstrained = 'PROJECTED GRADIENT';

%optimizer = 'IPOPT';

designVariable = 'MicroParams';
ub = 0.989;
lb = 0.011;
homegenizedVariablesComputer = 'ByVademecum';
% 
vademecumFileName = 'SuperEllipseQMax';
%vademecumFileName = 'SuperEllipseQ2';
%vademecumFileName = 'SuperEllipseQOptAnalytic';



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
% rho0 = 0.3;

line_search_initiator = 'INCREASING LAST STEP';
incrementFactor = 2.1;
%


kfrac = 2;
nsteps = 8;

plotting = true;
printing = true;
monitoring = true;
monitoring_interval = 2;
maxiter = 1600;


