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

%% DOMAIN DECOMPOSITION (1 SUBDOMAIN + DIRICHLET)
% Mesh building
mesh = data.mesh;

% Boundary mesh building
boundaryMesh = mesh.createBoundaryMesh();
boundaryDirichletMesh = boundaryMesh{1};

% Stiffness matrix computation
s.fValues = zeros(mesh.nnodes,mesh.ndim);
s.mesh    = mesh;
p1func = P1Function(s);

s.type     = 'ElasticStiffnessMatrix';
s.mesh     = mesh;
s.fun      = p1func;
s.material = data.material;
lhs = LHSintegrator.create(s);
K = lhs.compute();

% Mass matrix computation
trialDirichletMesh = P1Function.create(boundaryDirichletMesh.mesh,2);
testDirichletMesh = P1Function.create(boundaryDirichletMesh.mesh,2); 

s.type  = 'MassMatrix';
s.mesh  = boundaryDirichletMesh.mesh;
s.test  = testDirichletMesh;
s.trial = trialDirichletMesh;
lhs     = LHSintegrator.create(s);
C     = lhs.compute();

% DOFs definition
idxDirichletNodes = find(mesh.coord(:,1) == 0);
nDofDomainDirichlet = length(idxDirichletNodes)*mesh.ndim;
DofsDomainDirichlet = zeros(1,nDofDomainDirichlet);
for i=1:mesh.ndim
    DofsDomainDirichlet(i:2:(end-mesh.ndim+i)) = mesh.ndim*idxDirichletNodes+(-mesh.ndim+i);
end

idxNodes = 1:1:size(mesh.coord(),1);
nDofDomainMesh = length(idxNodes)*mesh.ndim;
DofsDomain = zeros(1,nDofDomainMesh);
for i=1:mesh.ndim
    DofsDomain(i:2:(end-mesh.ndim+i)) = mesh.ndim*idxNodes+(-mesh.ndim+i);
end
nDofTotal = max(DofsDomain);

if isa(testDirichletMesh,"P0Function") 
    nDofLagrange = boundaryDirichletMesh.mesh.nelem*testDirichletMesh.ndimf;
elseif isa(testDirichletMesh,"P1Function") 
    nDofLagrange = boundaryDirichletMesh.mesh.nnodes*testDirichletMesh.ndimf;
end
initialDofLagrange = max(nDofDomainMesh);
DofsLagrange = 1:1:nDofLagrange;
DofsLagrange = DofsLagrange + initialDofLagrange;
nDofTotal = nDofTotal + nDofLagrange;

% If lagrange multipliers use Dirac for their shape functions, use this:
%C = eye(nDofDomainDirichlet);

% Global LHS matrix assembly
GlobalLHS = sparse(nDofTotal,nDofTotal);
GlobalLHS(DofsDomain,DofsDomain) = K;
GlobalLHS(DofsDomainDirichlet,DofsLagrange) = C';
GlobalLHS(DofsLagrange,DofsDomainDirichlet) = C;

% Global RHS matrix assembly
GlobalRHS = sparse(nDofTotal,1);
F = zeros(nDofDomainMesh,1);
DofsNeumann = 2*data.bc.pointload(:,1)-2 + data.bc.pointload(:,2);
F(DofsNeumann) = data.bc.pointload(:,3);
GlobalRHS(1:nDofDomainMesh) = F;

% Solver (Direct)
Globalu = GlobalLHS\GlobalRHS;
Globalu = full(Globalu);

lambdaDirichlet = Globalu(nDofDomainMesh+1:end); 
u = Globalu(1:nDofDomainMesh);
u = [u(1:2:end-1), u(2:2:end)];

% Maximum error calculation
Error = max(max(abs(uTest.fValues-u)))

% Lagrange multipliers treatment
React = C*lambdaDirichlet;
React = [React(1:2:end-1), React(2:2:end)];
TotalReact = sum(React)

lambdaDirichlet = [lambdaDirichlet(1:2:end-1), lambdaDirichlet(2:2:end)];
%% PLOTS
% Displacement
s.fValues = uTest.fValues;
s.mesh    = mesh;
uTestPlot = P1Function(s);
uTestPlot.plot

s.fValues = u;
s.mesh    = mesh;
uPlot = P1Function(s);
uPlot.plot

% Lagrange multipliers (vectors)
if isa(testDirichletMesh,"P0Function") 
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
