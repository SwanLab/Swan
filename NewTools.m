
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

% Creating local domain 1 with unfitted mesh
sG.type            = 'Rectangle';
sG.xCoorCenter     = 0.5;
sG.yCoorCenter     = 0.5;
sG.xSide           = 1;
sG.ySide           = 1;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
uMesh1             = UnfittedMesh(sUm);
uMesh1.compute(lsFun.fValues);
figure;
plot(uMesh1);

% Creating local domain 2 with unfitted mesh
sG.type            = 'Rectangle';
sG.xCoorCenter     = 1.5;
sG.yCoorCenter     = 0.5;
sG.xSide           = 1;
sG.ySide           = 1;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
uMesh2             = UnfittedMesh(sUm);
uMesh2.compute(lsFun.fValues);
figure;
plot(uMesh2);

% Creating local perimeter constraint 1
sC.mesh       = mesh;
sC.uMesh      = uMesh1;
sC.epsilon    = 2*mesh.computeMeanCellSize();
sC.filter     = perFilter;
sC.value0     = 1;
sC.minEpsilon = mesh.computeMeanCellSize();
sC.target     = 1;
G2            = PerimeterConstraint(sC);
[J1,dJ1]      = G2.computeFunctionAndGradient(dens);

% Creating local perimeter constraint 2
sC.mesh       = mesh;
sC.uMesh      = uMesh2;
sC.epsilon    = 2*mesh.computeMeanCellSize();
sC.filter     = perFilter;
sC.value0     = 1.25;
sC.minEpsilon = mesh.computeMeanCellSize();
sC.target     = 1;
G2            = PerimeterConstraint(sC);
[J2,dJ2]      = G2.computeFunctionAndGradient(dens);

% Gradients sum
dJ = dJ1 + dJ2;
plot(dJ);