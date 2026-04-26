
clear;
clc;

mesh = TriangleMesh(1,1,50,50);
h    = mesh.computeMeanCellSize();

sG.type        = 'Circle';
sG.radius      = 0.35;
sG.xCoorCenter = 0.5;
sG.yCoorCenter = 0.5;
g              = GeometricalFunction(sG);
lsF            = g.computeLevelSetFunction(mesh);

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
uM                 = UnfittedMesh(sUm);
uM.compute(lsF.fValues);

chi = CharacteristicFunction.create(uM);

sF.mesh        = mesh;
sF.trial       = LagrangianFunction.create(mesh,1,'P1');
sF.senseVector = ConstantFunction.create([0,1],mesh);
sF.ovAngleDeg  = 45;
filter         = FilterOverhang(sF);
filter.updateEpsilon(4*h);

rhoEps = filter.compute(chi,3);