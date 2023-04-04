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
%filename = 'ArchTriFine';
%'Arch_quad_coarse';
%'Bridge_quad_coarse';
%'BridgeCool_Quadrilateral_Bilinear_Structured_Coarse';
%filename = 'Bridge';



ptype = 'MACRO';
initial_case = 'given';
m1 = 0.0101;
m2 = 0.0101;
cost = {'compliance'};
%cost = {'stressNorm'};
%cost = {'stressNorm','compliance'};
%cost = {'stressNorm','compliance'};
%weights = [1,10];
weights = 1;
constraint = {'volumeConstraint'};
filterType = 'PDE';
%filterType = 'P1';
constraint_case = 'EQUALITY';

Vfrac_initial = 0.3;
optimality_initial = 1e-5;
constr_initial = 1e-5;

Vfrac_final = 0.3;
optimality_final = 1e-5;
constr_final = 1e-5;

stressNormExponent_initial = 2;
stressNormExponent_final = 16;
% 
optimizer = 'DualNestedInPrimal';
%optimizer = 'AlternatingPrimalDual';
optimizerUnconstrained = 'PROJECTED GRADIENT';

%optimizer = 'MMA';
%optimizer = 'IPOPT';

designVariable = 'MicroParams';
ub = 0.989;
lb = 0.011;
homegenizedVariablesComputer = 'ByVademecum';
% 
%vademecumFileName = 'SuperEllipseQMax';
%vademecumFileName = 'SuperEllipseQ2';
vademecumFileName = 'SuperEllipseQOptAnalytic';



% designVariable = 'LevelSet';
% homegenizedVariablesComputer = 'ByInterpolation';
% method = 'SIMPALL';
% materialType = 'ISOTROPIC';
% initial_case = 'full';% optimizerUnconstrained = 'SLERP';



% designVariable = 'Density';
% homegenizedVariablesComputer = 'ByInterpolation';
% method = 'SIMPALL';
% materialType = 'ISOTROPIC';
% initial_case = 'full';
% rho0 = 0.3;

line_search_initiator = 'INCREASING LAST STEP';
incrementFactor = 1.05;
%


%kfrac = 2;
nsteps = 1;%17;

plotting = true;
printing = true;
monitoring = true;
monitoring_interval = 1;
maxiter = 800;

% 
% % % 
% isDirichletPartX = @(x) x > -1e-12 & x < 0.1;
% isDirichletPart1 = @(y) y > 0.20 & y < 0.30;
% isDirichletPart2 = @(y) y > 0.70 & y < 0.80;
% isDirichletPartY = @(y) isDirichletPart1(y) | isDirichletPart2(y);
% isDirichletPart = @(x,y) isDirichletPartX(x) & isDirichletPartY(y);
%  isNeumannPartX = @(x) x > (0.96) & x < (1+1e-12);
%  isNeumannPartY = @(y) y > 0.35 & y < 0.5;
%  
%  isCornerX = @(x) x > (0.4 - 0.05) & x < (0.4 +0.05);
%  isCornerY = @(y) y > (0.4 - 0.05) & y < (0.4 +0.05);
%  isCornerPart = @(x,y) isCornerX(x) & isCornerY(y);
% % 
%  isNeumannPart = @(x,y) isNeumannPartX(x) & isNeumannPartY(y);
% % iNotOptimizable = @(coord) isDirichletPart(coord(:,1),coord(:,2)) & ~isNeumannPart(coord(:,1),coord(:,2));
% 
% iNotOptimizable = @(coord) isCornerPart(coord(:,1),coord(:,2)) | isNeumannPart(coord(:,1),coord(:,2));
% 
% costDomainNotOptimizable       = {iNotOptimizable};
% constraintDomainNotOptimizable = {[]};
% 
% isDesignVariableFixed.nodes  = iNotOptimizable;
% isDesignVariableFixed.values = @(x) [m1*ones(size(x,1),1);m2*ones(size(x,1),1)];

