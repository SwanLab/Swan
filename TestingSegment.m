
clc;
clear;

m = TriangleMesh(1,1,100,100);

gPar.type        = 'Circle';
gPar.radius      = 0.3;
gPar.xCoorCenter = 0.5;
gPar.yCoorCenter = 0.5;
g                = GeometricalFunction(gPar);
phiFun           = g.computeLevelSetFunction(m);
phi              = phiFun.fValues;

sUm.backgroundMesh = m;
sUm.boundaryMesh   = m.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(phi);

funLS     = CharacteristicFunction.create(uMesh);

s.mesh  = m;
s.theta = 90;
s.alpha = 8;
s.beta  = 0;
filter  = NonLinearFilterSegment(s);
rhoEps  = filter.compute(funLS,3);