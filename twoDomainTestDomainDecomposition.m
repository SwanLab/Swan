clc,clear,close all

%% CONVENTIONAL SOLVER
% Create geometry
fileName = 'CantileverAxialLoad';
s.testName    = [fileName,'.m'];
s.problemCase = 'cantilever';
s.x1          = 1;
s.y1          = 0.5;
s.N           = 49;
s.M           = 100;
s.P           = 1;
s.DoF         = 2;

CantileverAxial = FEMInputWriter(s);
CantileverAxial.createTest();

% Create elasticity problem
a.fileName = fileName;
data = FemDataContainer(a);
fem = FEM.create(data);
fem.solve();
% fem.print(fileName);

% Get main results
uTest = fem.uFun(1,1);
e_test = fem.strainFun(1,1);
sig_test = fem.stressFun(1,1);

%% DOMAIN DECOMPOSITION (2 SUBDOMAINS)
% Mesh building
mesh = data.mesh;
xV = mesh.computeBaricenter();       
xBaricenter = xV(1,:);
yBaricenter = xV(2,:);
 
isLeft = xBaricenter <= max(mesh.coord(:,1))/2;
leftConnec = mesh.connec(isLeft,:);
s.connec = leftConnec;
s.coord = mesh.coord();
leftMesh = Mesh(s);
leftMesh = leftMesh.computeCanonicalMesh();

isRight = xBaricenter > max(mesh.coord(:,1))/2;
rightConnec = mesh.connec(isRight,:);
s.connec = rightConnec;
s.coord = mesh.coord;
rightMesh = Mesh(s);
rightMesh = rightMesh.computeCanonicalMesh();

% Boundary meshes building 
boundaryMesh = leftMesh.createBoundaryMesh();
boundaryDirichletMesh = boundaryMesh{1};
boundaryConnectionLeftMesh = boundaryMesh{2};
boundaryMesh = rightMesh.createBoundaryMesh();
boundaryConnectionRightMesh = boundaryMesh{1};
boundaryNeumannMesh = boundaryMesh{2};

% Stiffness matrices computation
% Left
s.fValues = zeros(leftMesh.nnodes, leftMesh.ndim);
s.mesh    = leftMesh;
p1left = P1Function(s);
s.type     = 'ElasticStiffnessMatrix';
s.mesh     = leftMesh;
s.fun      = p1left;
s.material = data.material;
lhs = LHSintegrator.create(s);
Kleft = lhs.compute();

% Right
s.fValues = zeros(rightMesh.nnodes, rightMesh.ndim);
s.mesh    = rightMesh;
p1right = P1Function(s);
s.type     = 'ElasticStiffnessMatrix';
s.mesh     = rightMesh;
s.fun      = p1right;
s.material = data.material;
lhs = LHSintegrator.create(s);
Kright = lhs.compute();

% Mass matrix computation
trialDirichletMesh = P1Function.create(boundaryDirichletMesh.mesh,2);
testDirichletMesh = P0Function.create(boundaryDirichletMesh.mesh,2);
s.type       = 'MassMatrix';
s.mesh       = boundaryDirichletMesh.mesh;
s.test       = testDirichletMesh;
s.trial      = trialDirichletMesh;
lhs          = LHSintegrator.create(s);
cDirichlet   = lhs.compute();

trialConnectionMesh = P1Function.create(boundaryConnectionLeftMesh.mesh,2);
testConnectionMesh = P0Function.create(boundaryConnectionLeftMesh.mesh,2);
s.type       = 'MassMatrix';
s.mesh       = boundaryConnectionLeftMesh.mesh;
s.test       = testDirichletMesh;
s.trial      = trialDirichletMesh;
lhs          = LHSintegrator.create(s);
cConnection  = lhs.compute();

% DOFs definition
nDofTotal = 0;
% Left
ndim = boundaryDirichletMesh.mesh.ndim;
idxDirichletNodes = find(ismember(leftMesh.coord,boundaryDirichletMesh.mesh.coord,'rows'));
nDofDirichlet = length(idxDirichletNodes)*ndim;
DofsDirichlet = zeros(1,nDofDirichlet);
for i=1:ndim
    DofsDirichlet(i:2:(end-ndim+i)) = ndim*idxDirichletNodes+(-ndim+i);
end

ndim = boundaryConnectionLeftMesh.mesh.ndim;
idxLeftConnectionNodes = find(ismember(leftMesh.coord,boundaryConnectionLeftMesh.mesh.coord,'rows'));
nDofConnection = length(idxLeftConnectionNodes)*ndim;
DofsLeftConnection = zeros(1,nDofConnection);
for i=1:ndim
    DofsLeftConnection(i:2:(end-ndim+i)) = ndim*idxLeftConnectionNodes+(-ndim+i);
end

ndim = leftMesh.ndim;
idxLeftNodes = 1:1:size(leftMesh.coord(),1);
nDofLeftDomain = length(idxLeftNodes)*ndim;
DofsLeftDomain = zeros(1,nDofLeftDomain);
for i=1:ndim
    DofsLeftDomain(i:2:(end-ndim+i)) = ndim*idxLeftNodes+(-ndim+i);
end
DofsLeftDomain = DofsLeftDomain + nDofTotal;
nDofTotal = max(DofsLeftDomain);

% Right Domain
ndim = boundaryNeumannMesh.mesh.ndim;
idxNeumannNodes = find(ismember(rightMesh.coord,boundaryNeumannMesh.mesh.coord,'rows'));
nDofNeumann = length(idxNeumannNodes)*boundaryNeumannMesh.mesh.ndim;
DofsNeumann = zeros(1,nDofNeumann);
for i=1:ndim
    DofsNeumann(i:2:(end-ndim+i)) = ndim*idxNeumannNodes+(-ndim+i);
end
DofsNeumann = DofsNeumann + nDofTotal;

ndim = boundaryConnectionLeftMesh.mesh.ndim;
idxRightConnectionNodes = find(ismember(rightMesh.coord,boundaryConnectionLeftMesh.mesh.coord,'rows'));
nDofConnection = length(idxRightConnectionNodes)*ndim;
DofsRightConnection = zeros(1,nDofConnection);
for i=1:ndim
    DofsRightConnection(i:2:(end-ndim+i)) = ndim*idxRightConnectionNodes+(-ndim+i);
end
DofsRightConnection = DofsRightConnection + nDofTotal;

ndim = rightMesh.ndim;
idxRightNodes = 1:1:size(rightMesh.coord(),1);
nDofRightDomain = length(idxRightNodes)*ndim;
DofsRightDomain = zeros(1,nDofRightDomain);
for i=1:ndim
    DofsRightDomain(i:2:(end-ndim+i)) = ndim*idxRightNodes+(-ndim+i);
end
DofsRightDomain = DofsRightDomain + nDofTotal;
nDofTotal = max(DofsRightDomain);

% Lagrange multipliers
if isa(testDirichletMesh,"P0Function") 
    nDofLagrangeDirichlet = boundaryDirichletMesh.mesh.nelem*testDirichletMesh.ndimf;
elseif isa(testDirichletMesh,"P1Function") 
    nDofLagrangeDirichlet = boundaryDirichletMesh.mesh.nnodes*testDirichletMesh.ndimf;
end
DofsLagrangeDirichlet = 1:1:nDofLagrangeDirichlet;
DofsLagrangeDirichlet = DofsLagrangeDirichlet + nDofTotal;
nDofTotal = max(DofsLagrangeDirichlet);

if isa(testConnectionMesh,"P0Function") 
    nDofLagrangeConnection = boundaryConnectionLeftMesh.mesh.nelem*testConnectionMesh.ndimf;
elseif isa(testConnectionMesh,"P1Function") 
    nDofLagrangeConnection = boundaryConnectionLeftMesh.mesh.nnodes*testConnectionMesh.ndimf;
end
DofsLagrangeConnection = 1:1:nDofLagrangeConnection;
DofsLagrangeConnection = DofsLagrangeConnection + nDofTotal;
nDofTotal = max(DofsLagrangeConnection);

% If lagrange multipliers use Dirac for their shape functions, use this:
%nDofLagrangeDirichlet = nDofDirichlet;
%nDofLagrangeConnection = nDofConnection;
%cDirichlet = eye(nDofDirichlet,nDofLagrangeDirichlet);
%cConnection = eye(nDofConnection, nDofLagrangeConnection);

% 7. Global LHS matrix assembly
GlobalLHS = sparse(nDofTotal, nDofTotal);
GlobalLHS(DofsLeftDomain,DofsLeftDomain) = Kleft;
GlobalLHS(DofsRightDomain,DofsRightDomain) = Kright;
GlobalLHS(DofsDirichlet,DofsLagrangeDirichlet) = cDirichlet';
GlobalLHS(DofsLeftConnection,DofsLagrangeConnection) = cConnection';
GlobalLHS(DofsRightConnection,DofsLagrangeConnection) = -cConnection';
GlobalLHS(DofsLagrangeDirichlet,DofsDirichlet) = cDirichlet;
GlobalLHS(DofsLagrangeConnection,DofsLeftConnection) = cConnection;
GlobalLHS(DofsLagrangeConnection,DofsRightConnection) = -cConnection;

% Global RHS matrix assembly
GlobalRHS = sparse(nDofTotal,1);
GlobalRHS(DofsNeumann) = data.bc.pointload(:,3);

% Solver 
Globalu = GlobalLHS\GlobalRHS;
Globalu = full(Globalu);

uLeft = Globalu(1:nDofLeftDomain);
uLeft = [uLeft(1:2:end-1), uLeft(2:2:end)];
uRight = Globalu(nDofLeftDomain+1:nDofLeftDomain + nDofRightDomain);
uRight = [uRight(1:2:end), uRight(2:2:end)];
lambdaDirichlet = Globalu(nDofLeftDomain+nDofRightDomain+1:nDofLeftDomain+nDofRightDomain+nDofLagrangeDirichlet);
lambdaConnection = Globalu(nDofLeftDomain+nDofRightDomain+nDofLagrangeDirichlet+1:nDofLeftDomain+nDofRightDomain+nDofLagrangeDirichlet+nDofLagrangeConnection);

% Maximum error calculation
ErrorConnectionDisp = max(max(abs(GlobalLHS(DofsLagrangeConnection,:)*Globalu)))

% Lagrange multipliers treatment
React = cDirichlet'*lambdaDirichlet;
React = [React(1:2:end-1), React(2:2:end)];
TotalReact = sum(React)

lambdaDirichlet = [lambdaDirichlet(1:2:end-1), lambdaDirichlet(2:2:end)];

%% PLOTS
% Subdomain meshes
figure
leftMesh.plot
figure
rightMesh.plot

% Displacement
s.fValues = uLeft;
s.mesh    = leftMesh;
uPlotLeft = P1Function(s);
uPlotLeft.plot

s.fValues = uRight;
s.mesh    = rightMesh;
uPlotRight = P1Function(s);
uPlotRight.plot

% Lagrange multipliers
if isa(testDirichletMesh,"P0Function") 
    % vectors
    x=(boundaryDirichletMesh.mesh.coord(:,1));
    x_avg= (x(1:end-1)+x(2:end))/2;
    
    y=(boundaryDirichletMesh.mesh.coord(:,2));
    y_avg= (y(1:end-1)+y(2:end))/2;
    
    figure
    quiver(x_avg,y_avg,lambdaDirichlet(:,1),zeros(size(lambdaDirichlet,1),1))
    figure
    quiver(x_avg,y_avg, zeros(size(lambdaDirichlet,1),1),lambdaDirichlet(:,2))
    figure
    quiver(x_avg,y_avg,lambdaDirichlet(:,1),lambdaDirichlet(:,2))

elseif isa(testDirichletMesh,"P1Function") 
    % function
    s.fValues = lambdaDirichlet(:,1);
    s.mesh    = boundaryDirichletMesh.mesh;
    LambdaDirichletXPlot = P1Function(s);
    LambdaDirichletXPlot.plot
    
    s.fValues = lambdaDirichlet(:,2);
    s.mesh    = boundaryDirichletMesh.mesh;
    LambdaDirichletYPlot = P1Function(s);
    LambdaDirichletYPlot.plot

    % vectors
    figure
    quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
           lambdaDirichlet(:,1),zeros(size(lambdaDirichlet,1),1))
    
    figure
    quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
           zeros(size(lambdaDirichlet,1),1),lambdaDirichlet(:,2))
    
    figure
    quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
           lambdaDirichlet(:,1),lambdaDirichlet(:,2))
end