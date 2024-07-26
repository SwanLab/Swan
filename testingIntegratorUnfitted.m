clear;
close all;

x1       = linspace(0,1,100);
x2       = linspace(0,1,100);
[xv,yv]  = meshgrid(x1,x2);
[F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
m.coord  = V(:,1:2);
m.connec = F;
mesh     = Mesh.create(m);
%mesh     = UnitQuadMesh(100,100);

gPar.type         = 'CircleInclusion';
gPar.radius       = 0.25;
gPar.xCoorCenter  = 0.5;
gPar.yCoorCenter  = 0.5;
g                 = GeometricalFunction(gPar);
phiFun            = g.computeLevelSetFunction(mesh);

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
uMesh              = UnfittedMesh(sUm);
uMesh.compute(phiFun.fValues);

chi = CharacteristicFunction.create(uMesh);

% % % Ex1 Scalar domain
ss.filterType   = 'PDE';
ss.mesh         = mesh;
ss.boundaryType = 'Neumann';
ss.metric       = 'Isotropy';
ss.trial = LagrangianFunction.create(mesh,1,'P1');
filter          = Filter.create(ss);

epsilon = 2*mesh.computeMeanCellSize(); %%
filter.updateEpsilon(epsilon); %%

rhoe          = filter.compute(chi,'QUADRATIC');

sUf.fun   = (1-rhoe);
sUf.uMesh = uMesh;
f         = UnfittedFunction(sUf);

int = Integrator.create('Unfitted',uMesh,'QUADRATIC');
P   = (2/epsilon)*int.compute(f);

% % % Ex2 Scalar boundary
chiB = CharacteristicFunction.createAtBoundary(uMesh);
totP = int.compute(chiB);

% % % Ex3 RHS domain
s.mesh     = uMesh;
s.type     = 'Unfitted';
s.quadType = 'QUADRATIC';
rhsInt     = RHSintegrator.create(s);
test       = LagrangianFunction.create(mesh, 1, 'P1');
intChiNi   = rhsInt.compute(chi,test);

% % % Ex4 RHS boundary
intChiBNi = rhsInt.compute(chiB,test);