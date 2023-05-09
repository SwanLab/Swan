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
s.N           = 20;
s.M           = 10;
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
u = fem.uFun{1,1};
% u.plot();
e = fem.strainFun{1,1};
% e.plot();
sig = fem.stressFun{1,1};
% sig.plot();


%% STEP 0: No domain decomposition (ONLY DIRICHLET)
% Mesh building
mesh = s.mesh;
bb.coord = mesh.coord(mesh.coord(:,1) == 0,:);
bb.connec = 1:9;
bb.connec = [bb.connec; bb.connec+1]';
bb.kFace = -1;
dirichletMesh =  Mesh(bb);

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

nDofLagrange = nDofDomainDirichlet;
initialDofLagrange = max(DofsDomain)+1;
DofsLagrange = initialDofLagrange:1:(initialDofLagrange+nDofDomainDirichlet-1);

C = eye(nDofDomainDirichlet);

% Global LHS matrix assembly
nDofTotal = nDofLagrange + nDofDomainMesh;
GlobalLHS = sparse(nDofTotal);

GlobalLHS(DofsDomain,DofsDomain) = K;
GlobalLHS(DofsDomainDirichlet,DofsLagrange) = C;
GlobalLHS(DofsLagrange,DofsDomainDirichlet) = C';

% Global RHS matrix assembly
GlobalRHS = zeros(nDofTotal,1);
F = zeros(nDofDomainMesh,1);
DofsNeumann = 2*s.bc.pointload(:,1)-2 + s.bc.pointload(:,2);
F(DofsNeumann) = s.bc.pointload(:,3);
GlobalRHS(1:nDofDomainMesh) = F;

% Solver (Direct)
u_dd = GlobalLHS\GlobalRHS;
u_dd = u_dd(1:nDofDomainMesh);
uDD = [u_dd(1:2:end-1), u_dd(2:2:end)];
Error = abs(u.fValues-uDD);
max(max(Error))

zz.fValues = u.fValues;
zz.mesh    = mesh;
plotUFD = P1Function(zz);

zz.fValues = uDD;
zz.mesh    = mesh;
plotUDD = P1Function(zz);

%plots
plotUFD.plot
plotUDD.plot
Error = max(max(abs(u.fValues-uDD)))

%% STEP 1.1: Mesh decomposition
mesh = s.mesh;
xV = mesh.computeBaricenter();       
x_baricenter = xV(1,:);
y_baricenter = xV(2,:);
 
isLeft = x_baricenter <= max(mesh.coord(:,1))/2;
leftConnec = mesh.connec(isLeft,:);
isRight = x_baricenter > max(mesh.coord(:,1))/2;
rightConnec = mesh.connec(isRight,:);
 
ss.connec = leftConnec;
ss.coord = mesh.coord;
leftMesh = Mesh(ss);
leftMesh = leftMesh.computeCanonicalMesh();
boundaryLeftMesh = leftMesh.createBoundaryMesh();

ss.connec = rightConnec;
ss.coord = mesh.coord;
rightMesh = Mesh(ss);
rightMesh = rightMesh.computeCanonicalMesh();
boundaryRightMesh = rightMesh.createBoundaryMesh();
 
% % Boundary mesh 
% bb.coord = leftMesh.coord(leftMesh.coord(:,1) == 0.5,:);
% bb.connec = [1:9];
% bb.connec = [bb.connec; bb.connec+1]';
% bb.kFace = -1;
% leftBoundaryMesh = Mesh(bb);
% 
% bb.coord = rightMesh.coord(rightMesh.coord(:,1) == 0.5,:);
% bb.connec = [1:9];
% bb.connec = [bb.connec; bb.connec+1]';
% bb.kFace = -1;
% rightBoundaryMesh = Mesh(bb);

 
% Plots
figure
leftMesh.plot
figure
rightMesh.plot

fValuesLeft = u.fValues(unique(leftMesh.connec()),:);
fValuesRight = u.fValues(unique(rightMesh.connec()),:);

zz.fValues = fValuesLeft;
zz.mesh    = leftMesh;
p1left = P1Function(zz);

zz.fValues = fValuesRight;
zz.mesh    = rightMesh;
p1right = P1Function(zz);

% "Mass matrix definition (P1 for displacement & P0 for lagrange multipliers)
% sss.mesh    = leftBoundaryMesh;
% nnode       = size(leftBoundaryMesh.coord,1);
% ndof        = 1;
% sss.fValues = ones(nnode*ndof,1);
% P1LeftBMesh = P1Function(sss);
% P1LeftBMesh.plot
% 
% sss.mesh    = leftBoundaryMesh;
% nnode       = size(leftBoundaryMesh.coord,1);
% ndof        = 1;
% nelem = nnode - 1; 
% sss.fValues = ones(nelem*ndof,1);
% P0LeftBMesh = P0Function(sss);
% P0LeftBMesh.plot

% ... inside LHS
% quad        = Quadrature.create(leftBoundaryMesh,'LINEAR');
% N           = P1LeftBMesh.computeShapeFunctions(quad);
%

%% STEP 1.2: Stiffness matrix
% % ...
% sLHS.type     = 'ElasticStiffnessMatrix';
% sLHS.mesh     = leftMesh;
% sLHS.fun      = p1left;
% sLHS.material = s.material;
% lhs = LHSintegrator.create(sLHS);
% Kleft = lhs.compute();
% 
% sLHS.type     = 'ElasticStiffnessMatrix';
% sLHS.mesh     = rightMesh;
% sLHS.fun      = p1right;
% sLHS.material = s.material;
% lhs = LHSintegrator.create(sLHS);
% Kright = lhs.compute();

%% STEP 4: "Mass matrices" definition and assembly







