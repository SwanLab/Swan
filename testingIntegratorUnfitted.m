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

chi = CharacteristicFunction.create(uMesh);