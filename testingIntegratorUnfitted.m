clear;
close all;

x1       = linspace(0,1,100);
x2       = linspace(0,1,100);
[xv,yv]  = meshgrid(x1,x2);
[F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
m.coord  = V(:,1:2);
m.connec = F;
mesh     = Mesh.create(m);

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

chi           = CharacteristicFunction.create(uMesh);
ss.filterType   = 'PDE';
ss.mesh         = mesh;
ss.boundaryType = 'Neumann';
ss.metric       = 'Isotropy';
ss.trial = LagrangianFunction.create(mesh,1,'P1');
filter          = Filter.create(ss);

epsilon = 2*mesh.computeMeanCellSize(); %%
filter.updateEpsilon(epsilon); %%

rhoe          = filter.compute(chi,'QUADRATICMASS');

sUf.fun   = (1-rhoe);
sUf.uMesh = uMesh;
f         = UnfittedFunction(sUf);
uF        = f.unfittedMeshFunction;

int = Integrator.create('Unfitted',uMesh,'QUADRATICMASS');
P   = (2/epsilon)*int.compute(uF);