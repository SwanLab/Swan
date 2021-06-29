filename = 'CantileverSquareSym';%'CantileverSquareNewFine';%'LshapeTriFine';%'CantileverSquareNew';'ArchTriFine';%'CantileverSquare';%'ArchTriFine';%'CantileverSquare';%'Lshape';%'LshapeTriFine';%'CantileverSquareSmall';%'';%'Lshape';'CantileverSquare';'Lshape';%'CantileverSquare';%'LshapeFine';%'Bridge_quad_coarse';%'Arch_quad_coarse';%'BridgeCool_Quadrilateral_Bilinear_Structured_Coarse';%'Bridge';%'CantileverSquareSmall';
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
%vademecumFileName = 'SuperEllipseQMax';
%vademecumFileName = 'SuperEllipseQ2';
vademecumFileName = 'SuperEllipseQOptAnalytic';
% 
% designVariable = 'Density';
% homegenizedVariablesComputer = 'ByInterpolation';
% method = 'SIMPALL';
% materialType = 'ISOTROPIC';

line_search_initiator = 'INCREASING LAST STEP';
incrementFactor = 1.95;
%

%kfrac = 2;
nsteps = 8;%17;

plotting = true;
printing = true;
monitoring = true;
monitoring_interval = 3;
maxiter = 800;



% % % 
isDirichletPartX = @(x) x > -1e-12 & x < 0.1;
%isDirichletPart1 = @(y) y > 0.20 & y < 0.30;
isDirichletPart1 = @(y) y > 0.30 & y < 0.50;
isDirichletPartY = @(y) isDirichletPart1(y) ;%| isDirichletPart2(y);
isDirichletPart = @(x,y) isDirichletPartX(x) & isDirichletPartY(y);
isNeumannPartX = @(x) x > (2-0.1) & x < (2+1e-12);
isNeumannPartY = @(y) y > 0.00 & y < 0.1;
isNeumannPart = @(x,y) isNeumannPartX(x) & isNeumannPartY(y);

iNotOptimizable = @(coord) isDirichletPart(coord(:,1),coord(:,2)) | isNeumannPart(coord(:,1),coord(:,2));

costDomainNotOptimizable       = {iNotOptimizable};
constraintDomainNotOptimizable = {[]};

isDesignVariableFixed.nodes  = iNotOptimizable;
isDesignVariableFixed.values = @(x) [m1*ones(size(x,1),1);m2*ones(size(x,1),1)];
