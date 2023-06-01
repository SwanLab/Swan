% Testing domain decomposition for reduced order modeling

clear
clc
close all

%% STEP 1: Compute full domain

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
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();
% fem.print(fileName);

% Get main results
u = fem.uFun(1,1);
% u.plot();
e = fem.strainFun(1,1);
% e.plot();
sig = fem.stressFun(1,1);
% sig.plot();


%% STEP 0: No domain decomposition (ONLY DIRICHLET)
% Mesh building
mesh = s.mesh;
boundaryMesh = mesh.createBoundaryMesh();
boundaryDirichletMesh = boundaryMesh{1};

% Domain's stiffness matrix
zz.fValues = u.fValues;
zz.mesh    = mesh;
p1func = P1Function(zz);

sLHS.type     = 'ElasticStiffnessMatrix';
sLHS.mesh     = mesh;
sLHS.fun      = p1func;
sLHS.material = s.material;
lhs = LHSintegrator.create(sLHS);
K = lhs.compute();

% 'Mass' matrices assembly
trialDirichletMesh = P1Function.create(boundaryDirichletMesh.mesh,2);
testDirichletMesh = P1Function.create(boundaryDirichletMesh.mesh,2);

ss.type  = 'MassMatrix';
ss.mesh  = boundaryDirichletMesh.mesh;
ss.test  = testDirichletMesh;
ss.trial = trialDirichletMesh;
lhs     = LHSintegrator.create(ss);
C     = lhs.compute();

idxNodes = 1:1:size(mesh.coord(),1);
nDofDomainMesh = length(idxNodes)*2;
DofsDomain = zeros(1,nDofDomainMesh);
DofsDomain(1:2:end-1) = 2*idxNodes-1;
DofsDomain(2:2:end) = 2*idxNodes;

idxDirichletNodes = find(mesh.coord(:,1) == 0);
nDofDomainDirichlet = length(idxDirichletNodes)*2;
DofsDomainDirichlet = zeros(1,nDofDomainDirichlet);
DofsDomainDirichlet(1:2:end-1) = 2*idxDirichletNodes-1;
DofsDomainDirichlet(2:2:end) = 2*idxDirichletNodes;

nDofLagrange = boundaryDirichletMesh.mesh.nnodes*testDirichletMesh.ndimf;
initialDofLagrange = max(nDofDomainMesh)+1;
DofsLagrange = initialDofLagrange:1:(initialDofLagrange+nDofLagrange-1);

%C = eye(nDofDomainDirichlet);

% Global LHS matrix assembly
nDofTotal = nDofLagrange + nDofDomainMesh;
GlobalLHS = sparse(nDofTotal);

GlobalLHS(DofsDomain,DofsDomain) = K;
GlobalLHS(DofsDomainDirichlet,DofsLagrange) = C';
GlobalLHS(DofsLagrange,DofsDomainDirichlet) = C;

% Global RHS matrix assembly
GlobalRHS = zeros(nDofTotal,1);
F = zeros(nDofDomainMesh,1);
DofsNeumann = 2*s.bc.pointload(:,1)-2 + s.bc.pointload(:,2);
F(DofsNeumann) = s.bc.pointload(:,3);
GlobalRHS(1:nDofDomainMesh) = F;

% Solver (Direct)
u_dd = GlobalLHS\GlobalRHS;
u_disp = u_dd(1:nDofDomainMesh);
uDD = [u_disp(1:2:end-1), u_disp(2:2:end)];

Error = max(max(abs(u.fValues-uDD)))

lambdaDirichlet = u_dd(nDofDomainMesh+1:end);
React = C*lambdaDirichlet
sum(React)
lambdaDirichlet = [React(1:2:end-1), React(2:2:end)];

%plots
zz.fValues = u.fValues;
zz.mesh    = mesh;
plotUFD = P1Function(zz);
plotUFD.plot

zz.fValues = uDD;
zz.mesh    = mesh;
plotUDD = P1Function(zz);
plotUDD.plot

%plot Lagrange values
zz.fValues = lambdaDirichlet(:,1);
zz.mesh    = boundaryDirichletMesh.mesh;
plotX = P1Function(zz);
plotX.plot

zz.fValues = lambdaDirichlet(:,2);
zz.mesh    = boundaryDirichletMesh.mesh;
plotY = P1Function(zz);
plotY.plot

%%% P1 Plot
figure
quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
       lambdaDirichlet(:,1),zeros(size(lambdaDirichlet,1),1))
figure
quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
       zeros(size(lambdaDirichlet,1),1),lambdaDirichlet(:,2))
figure
quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
       lambdaDirichlet(:,1),lambdaDirichlet(:,2))

%%% P0 plot
% x=(boundaryDirichletMesh.mesh.coord(:,1));
% x_avg= (x(1:end-1)+x(2:end))/2;
% 
% y=(boundaryDirichletMesh.mesh.coord(:,2));
% y_avg= (y(1:end-1)+y(2:end))/2;
% 
% figure
% quiver(x_avg,y_avg,lambdaDirichlet(:,1),zeros(size(lambdaDirichlet,1),1))
% figure
% quiver(x_avg,y_avg, zeros(size(lambdaDirichlet,1),1),lambdaDirichlet(:,2))
% figure
% quiver(x_avg,y_avg,lambdaDirichlet(:,1),lambdaDirichlet(:,2))

%% STEP 1.1: Mesh decomposition
% 1. Create subdomains
mesh = s.mesh;
xV = mesh.computeBaricenter();       
x_baricenter = xV(1,:);
y_baricenter = xV(2,:);
 
isLeft = x_baricenter <= max(mesh.coord(:,1))/2;
leftConnec = mesh.connec(isLeft,:);
isRight = x_baricenter > max(mesh.coord(:,1))/2;
rightConnec = mesh.connec(isRight,:);
 
% 2. Create subdomains' meshes
ss.connec = leftConnec;
ss.coord = mesh.coord;
leftMesh = Mesh(ss);
leftMesh = leftMesh.computeCanonicalMesh();

ss.connec = rightConnec;
ss.coord = mesh.coord;
rightMesh = Mesh(ss);
rightMesh = rightMesh.computeCanonicalMesh();

% Plots
figure
leftMesh.plot
figure
rightMesh.plot

% 3. Create boundaries' meshes 
boundaryMesh = leftMesh.createBoundaryMesh();
boundaryDirichletMesh = boundaryMesh{1};
boundaryConnectionLeftMesh = boundaryMesh{2};
boundaryMesh = rightMesh.createBoundaryMesh();
boundaryConnectionRightMesh = boundaryMesh{1};

% Plots
figure
boundaryDirichletMesh.mesh.plot
figure
boundaryConnectionLeftMesh.mesh.plot

% 4. Compute stiffness matrices
fValuesLeft = zeros(leftMesh.nnodes, leftMesh.ndim);
fValuesRight = zeros(rightMesh.nnodes, rightMesh.ndim);

zz.fValues = fValuesLeft;
zz.mesh    = leftMesh;
p1left = P1Function(zz);

zz.fValues = fValuesRight;
zz.mesh    = rightMesh;
p1right = P1Function(zz);

sLHS.type     = 'ElasticStiffnessMatrix';
sLHS.mesh     = leftMesh;
sLHS.fun      = p1left;
sLHS.material = s.material;
lhs = LHSintegrator.create(sLHS);
Kleft = lhs.compute();

sLHS.type     = 'ElasticStiffnessMatrix';
sLHS.mesh     = rightMesh;
sLHS.fun      = p1right;
sLHS.material = s.material;
lhs = LHSintegrator.create(sLHS);
Kright = lhs.compute();

% 5.1. Dirichlet matrix definition [Cd] (P1 for displacement & P0 for lagrange multipliers)
% ss.mesh      = boundaryDirichletMesh.mesh;
% nnode       = boundaryDirichletMesh.mesh.nnodes;
% ndim        = boundaryDirichletMesh.mesh.ndim;
% ss.fValues   = zeros(nnode,ndim);
P1DirichletMesh = P1Function.create(boundaryDirichletMesh.mesh,2);

% ss.mesh      = boundaryDirichletMesh.mesh;
% nelem       = boundaryDirichletMesh.mesh.nelem;
% ndim        = boundaryDirichletMesh.mesh.ndim;
% ss.fValues   = zeros(nelem,ndim);
P0DirichletMesh = P0Function.create(boundaryDirichletMesh.mesh,2);

ss.type  = 'MassMatrix';
ss.mesh  = boundaryDirichletMesh.mesh;
ss.test  = P0DirichletMesh;
ss.trial = P1DirichletMesh;
lhs     = LHSintegrator.create(ss);
cDirichlet     = lhs.compute();


% 6. DOFs allocation
% Left Domain
idxDirichletNodes = find(ismember(leftMesh.coord,boundaryDirichletMesh.mesh.coord,'rows'));
nDofDirichlet = length(idxDirichletNodes)*2;
DofsDirichlet = zeros(1,nDofDirichlet);
DofsDirichlet(1:2:end-1) = 2*idxDirichletNodes-1;
DofsDirichlet(2:2:end) = 2*idxDirichletNodes;

idxLeftConnectionNodes = find(ismember(leftMesh.coord,boundaryConnectionLeftMesh.mesh.coord,'rows'));
nDofConnection = length(idxLeftConnectionNodes)*2;
DofsLeftConnection = zeros(1,nDofConnection);
DofsLeftConnection(1:2:end-1) = 2*idxLeftConnectionNodes-1;
DofsLeftConnection(2:2:end) = 2*idxLeftConnectionNodes;

idxLeftNodes = 1:1:size(leftMesh.coord(),1);
nDofLeftDomain = length(idxLeftNodes)*2;
DofsLeftDomain = zeros(1,nDofLeftDomain);
DofsLeftDomain(1:2:end-1) = 2*idxLeftNodes - 1;
DofsLeftDomain(2:2:end) = 2*idxLeftNodes;

% Right Domain
idxNeumannNodes = find(rightMesh.coord(:,1) == 1);
nDofNeumann = length(idxNeumannNodes)*2;
DofsNeumann = zeros(1,nDofNeumann);
DofsNeumann(1:2:end-1) = 2*idxNeumannNodes-1;
DofsNeumann(2:2:end) = 2*idxNeumannNodes;
DofsNeumann = DofsNeumann + max(DofsLeftDomain);

idxRightConnectionNodes = find(ismember(rightMesh.coord,boundaryConnectionLeftMesh.mesh.coord,'rows'));
DofsRightConnection = zeros(1,nDofConnection);
DofsRightConnection(1:2:end-1) = 2*idxRightConnectionNodes-1;
DofsRightConnection(2:2:end) = 2*idxRightConnectionNodes;
DofsRightConnection = DofsRightConnection + max(DofsLeftDomain);

idxRightNodes = 1:1:size(rightMesh.coord(),1);
nDofRightDomain = length(idxRightNodes)*2;
DofsRightDomain = zeros(1,nDofRightDomain);
DofsRightDomain(1:2:end-1) = 2*idxRightNodes - 1;
DofsRightDomain(2:2:end) = 2*idxRightNodes;
DofsRightDomain = DofsRightDomain + max(DofsLeftDomain);

% Lagrange multipliers' mesh
nDofLagrangeDirichlet = boundaryDirichletMesh.mesh.nelem*P0DirichletMesh.ndimf;
initialDofLagrange = max(DofsRightDomain)+1;
DofsLagrangeDirichlet = initialDofLagrange:1:(initialDofLagrange+nDofLagrangeDirichlet-1);

nDofLagrangeConnection = nDofConnection; % nDofRightConnection must be the same
initialDofLagrange = max(DofsLagrangeDirichlet)+1;
DofsLagrangeConnection = initialDofLagrange:1:(initialDofLagrange+nDofConnection-1);

nDofTotal = max(DofsLagrangeConnection);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cDirichlet = eye(nDofDirichlet,nDofLagrangeDirichlet);
cConnection = eye(nDofConnection, nDofLagrangeConnection);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 7. Global Assembly
GlobalLHS = sparse(nDofTotal, nDofTotal);
GlobalLHS(DofsLeftDomain,DofsLeftDomain) = Kleft;
GlobalLHS(DofsRightDomain,DofsRightDomain) = Kright;
GlobalLHS(DofsDirichlet,DofsLagrangeDirichlet) = cDirichlet';
GlobalLHS(DofsLeftConnection,DofsLagrangeConnection) = cConnection;
GlobalLHS(DofsRightConnection,DofsLagrangeConnection) = -cConnection;
GlobalLHS(DofsLagrangeDirichlet,DofsDirichlet) = cDirichlet;
GlobalLHS(DofsLagrangeConnection,DofsLeftConnection) = cConnection;
GlobalLHS(DofsLagrangeConnection,DofsRightConnection) = -cConnection;

GlobalRHS = zeros(nDofTotal,1);
GlobalRHS(DofsNeumann) = s.bc.pointload(:,3);



% 8. Solver 
u = GlobalLHS\GlobalRHS;
uLeft = u(1:nDofLeftDomain);
uLeft = [uLeft(1:2:end-1), uLeft(2:2:end)];

uRight = u(nDofLeftDomain+1:nDofLeftDomain + nDofRightDomain);
uRight = [uRight(1:2:end), uRight(2:2:end)];

lambdaDirichlet = u(nDofLeftDomain+nDofRightDomain+1:nDofLeftDomain+nDofRightDomain+nDofLagrangeDirichlet);
lambdaDirichlet = [lambdaDirichlet(1:2:end-1), lambdaDirichlet(2:2:end)];

%plot (subdomains)
% zz.fValues = uLeft;
% zz.mesh    = leftMesh;
% plotLeft = P1Function(zz);
% plotLeft.plot
% 
% zz.fValues = uRight;
% zz.mesh    = rightMesh;
% plotRight = P1Function(zz);
% plotRight.plot

%plot Lagrange values
zz.fValues = lambdaDirichlet(:,1);
zz.mesh    = boundaryDirichletMesh.mesh;
plotX = P1Function(zz);
plotX.plot

zz.fValues = lambdaDirichlet(:,2);
zz.mesh    = boundaryDirichletMesh.mesh;
plotY = P1Function(zz);
plotY.plot

figure
quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
       lambdaDirichlet(:,1),zeros(size(lambdaDirichlet,1),1))

figure
quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
       zeros(size(lambdaDirichlet,1),1),lambdaDirichlet(:,2))

figure
quiver(boundaryDirichletMesh.mesh.coord(:,1),boundaryDirichletMesh.mesh.coord(:,2),...
       lambdaDirichlet(:,1),lambdaDirichlet(:,2))

% plot (combined)
%uGlobal = [uLeft; uRight];

%test
Error = max(max(abs(GlobalLHS(DofsLagrangeConnection,:)*u)))



