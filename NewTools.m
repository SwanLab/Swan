
clear;
close all;
clc;

% Domain example
mesh = TriangleMesh(2,1,80,40);

% Design variable
sG.type        = 'SmoothRectangle';
sG.xCoorCenter = 1;
sG.yCoorCenter = 0.5;
sG.xSide       = 1.25;
sG.ySide       = 0.5;
sG.pnorm       = 8;
g              = GeometricalFunction(sG);
lsFun          = g.computeLevelSetFunction(mesh);
sf.fValues     = 1-heaviside(lsFun.fValues);
sf.mesh        = mesh;
sf.order       = 'P1';
s.fun          = LagrangianFunction(sf);
s.mesh         = mesh;
s.type         = 'Density';
s.plotting     = false;
dens           = DesignVariable.create(s);
plot(dens.fun);

% Filter PDE for the perimeter
sfi.trial   = LagrangianFunction.create(mesh,1,'P1');
sfi.mesh    = mesh;
sfi.LHStype = 'StiffnessMass';
perFilter   = FilterPDE(sfi);

% Creating constraint of interest
sF.mesh               = mesh;
sF.epsilon            = 2*mesh.computeMeanCellSize();
sF.minEpsilon         = 2*mesh.computeMeanCellSize();
sF.perimeterTargetAbs = 1;
sF.filter             = perFilter;
globConstr            = PerimeterConstraint(sF);

% Creating local domain with unfitted mesh
sG.type            = 'Rectangle';
sG.xCoorCenter     = 0.5;
sG.yCoorCenter     = 0.5;
sG.xSide           = 1;
sG.ySide           = 1;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
uMesh              = UnfittedMesh(sUm);
uMesh.compute(lsFun.fValues);
figure;
plot(uMesh);

% Creating local constraint functional
sC.mesh         = mesh;
sC.constraint   = globConstr;
sC.unfittedMesh = uMesh;
sC.epsilon      = 2*mesh.computeMeanCellSize();
sC.target       = 1;
sC.filterPer    = perFilter;
sC.value0       = 1;
G               = LocalConstraint(sC);
[J,dJ]          = G.computeFunctionAndGradient(dens);
plot(dJ);