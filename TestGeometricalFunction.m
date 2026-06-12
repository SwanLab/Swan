clc,clear,close all

file = 'Tetradecahedron';
a.fileName = file;
fem = FemDataContainer(a);
%baseMesh = TetraMesh(1, 1, 1, 20, 20, 20);
baseMesh = fem.mesh;

gPar.type = 'Tetradecahedron'; %'Octahedron';
gPar.xCoorCenter = 0.5;
gPar.yCoorCenter = 0.5;
gPar.zCoorCenter = 0.5;
gPar.radius = 1e-5;
g         = GeometricalFunction(gPar);
phiFun    = g.computeLevelSetFunction(baseMesh);
lsCircle  = phiFun.fValues;
ls = -lsCircle;


sUm.backgroundMesh = baseMesh;
sUm.boundaryMesh   = baseMesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(ls);
uMesh.plot

ls = CharacteristicFunction.create(uMesh);
s.trial = LagrangianFunction.create(baseMesh,1,'P1');
s.mesh = baseMesh;
f = FilterLump(s);
lsf = f.compute(ls,2);
