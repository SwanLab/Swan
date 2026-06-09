close all
baseMesh = TetraMesh(1, 1, 1, 20, 20, 20);

gPar.type = 'Octahedron';%'Tetradecahedron';
gPar.xCoorCenter = 0.5;
gPar.yCoorCenter = 0.5;
gPar.zCoorCenter = 0.5;
gPar.radius = 3/4;
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
