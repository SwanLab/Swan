
clear;
clc;

mesh = TriangleMesh(1,1,100,100);
h    = mesh.computeMeanCellSize();

s.type        = 'Rectangle';
s.xSide       = 0.75;
s.ySide       = 0.25;
s.xCoorCenter = 0.5;
s.yCoorCenter = 0.75;
g             = GeometricalFunction(s);
ls1           = g.computeLevelSetFunction(mesh);

s.type        = 'Rectangle';
s.xSide       = 0.25;
s.ySide       = 0.75;
s.xCoorCenter = 0.5;
s.yCoorCenter = 0.5;
g             = GeometricalFunction(s);
ls2           = g.computeLevelSetFunction(mesh);

lsVal = min(ls1.fValues,ls2.fValues);

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
uMesh = UnfittedMesh(sUm);
uMesh.compute(lsVal);

chi = CharacteristicFunction.create(uMesh);

sF.mesh = mesh;
sF.trial      = LagrangianFunction.create(mesh,1,'P1');
fLump = FilterLump(sF);
xLump = fLump.compute(chi,2);

sF.filterType = 'PDE';
fPDE          = Filter.create(sF);
fPDE.updateEpsilon(3*h);
xPDE          = fPDE.compute(chi,2);

CAnisotropic = [tand(85), 0; 0, 1/tand(85)];
aniAlphaDeg = 90;
R = [cosd(aniAlphaDeg),-sind(aniAlphaDeg)
    sind(aniAlphaDeg), cosd(aniAlphaDeg)];
CGlobal = R*CAnisotropic*R';
sF.boundaryType = 'Neumann';
sF.metric       = 'Anisotropy';
sF.A            = ConstantFunction.create(CGlobal,mesh);
fAni            = Filter.create(sF);
fAni.updateEpsilon(3*h);
xAni = fAni.compute(chi,2);

sF.alpha = 4;
sF.beta  = 0;
sF.theta = 90;
sF.tol0  = 1e-6;
fSeg     = NonLinearFilterSegment(sF);
fSeg.updateEpsilon(3*h);
xSeg = fSeg.compute(chi,2);