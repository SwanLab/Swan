%% Pure Perimeter

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


%% Min compliance MBB

clear;
close all;
clc;

% Domain and design variable
load('NumericalExperiments/MBB/DesVarDensityMinCompliance.mat','d');
mesh = d.fun.mesh;

% Filter PDE for the perimeter
sfi.trial   = LagrangianFunction.create(mesh,1,'P1');
sfi.mesh    = mesh;
sfi.LHStype = 'StiffnessMass';
perFilter   = FilterPDE(sfi);

% Creating global domain
sG.type            = 'Full';
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
Omega              = UnfittedMesh(sUm);
Omega.compute(lsFun.fValues);

% Creating local domain Left
sG.type            = 'Rectangle';
sG.xCoorCenter     = 1.5;
sG.yCoorCenter     = 0.5;
sG.xSide           = 3;
sG.ySide           = 1;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaL             = UnfittedMesh(sUm);
OmegaL.compute(lsFun.fValues);

% Creating local domain Right
sG.type            = 'Rectangle';
sG.xCoorCenter     = 4.5;
sG.yCoorCenter     = 0.5;
sG.xSide           = 3;
sG.ySide           = 1;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaR             = UnfittedMesh(sUm);
OmegaR.compute(lsFun.fValues);

% Creating local domain LeftUp
sG.type            = 'Rectangle';
sG.xCoorCenter     = 1.5;
sG.yCoorCenter     = 0.75;
sG.xSide           = 3;
sG.ySide           = 0.5;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaLU            = UnfittedMesh(sUm);
OmegaLU.compute(lsFun.fValues);

% Creating local domain LeftDown
sG.type            = 'Rectangle';
sG.xCoorCenter     = 1.5;
sG.yCoorCenter     = 0.25;
sG.xSide           = 3;
sG.ySide           = 0.5;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaLD            = UnfittedMesh(sUm);
OmegaLD.compute(lsFun.fValues);

% Creating local domain RightUp
sG.type            = 'Rectangle';
sG.xCoorCenter     = 4.5;
sG.yCoorCenter     = 0.75;
sG.xSide           = 3;
sG.ySide           = 0.5;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaRU            = UnfittedMesh(sUm);
OmegaRU.compute(lsFun.fValues);

% Creating local domain RightDown
sG.type            = 'Rectangle';
sG.xCoorCenter     = 4.5;
sG.yCoorCenter     = 0.25;
sG.xSide           = 3;
sG.ySide           = 0.5;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaRD            = UnfittedMesh(sUm);
OmegaRD.compute(lsFun.fValues);

% Creating omega perimeter
sP.mesh       = mesh;
sP.uMesh      = Omega;
sP.epsilon    = mesh.computeMeanCellSize();
sP.filter     = perFilter;
sP.value0     = 1;
P             = PerimeterFunctional(sP);
[POmega,~]    = P.computeFunctionAndGradient(d);

% Creating omegaL perimeter
sP.mesh       = mesh;
sP.uMesh      = OmegaL;
sP.epsilon    = mesh.computeMeanCellSize();
sP.filter     = perFilter;
sP.value0     = 1;
P             = PerimeterFunctional(sP);
[POmegaL,~]    = P.computeFunctionAndGradient(d);

% Creating omegaR perimeter
sP.mesh       = mesh;
sP.uMesh      = OmegaR;
sP.epsilon    = mesh.computeMeanCellSize();
sP.filter     = perFilter;
sP.value0     = 1;
P             = PerimeterFunctional(sP);
[POmegaR,~]    = P.computeFunctionAndGradient(d);

% Creating omegaLU perimeter
sP.mesh       = mesh;
sP.uMesh      = OmegaLU;
sP.epsilon    = mesh.computeMeanCellSize();
sP.filter     = perFilter;
sP.value0     = 1;
P             = PerimeterFunctional(sP);
[POmegaLU,~]    = P.computeFunctionAndGradient(d);

% Creating omegaLD perimeter
sP.mesh       = mesh;
sP.uMesh      = OmegaLD;
sP.epsilon    = mesh.computeMeanCellSize();
sP.filter     = perFilter;
sP.value0     = 1;
P             = PerimeterFunctional(sP);
[POmegaLD,~]    = P.computeFunctionAndGradient(d);

% Creating omegaRU perimeter
sP.mesh       = mesh;
sP.uMesh      = OmegaRU;
sP.epsilon    = mesh.computeMeanCellSize();
sP.filter     = perFilter;
sP.value0     = 1;
P             = PerimeterFunctional(sP);
[POmegaRU,~]    = P.computeFunctionAndGradient(d);

% Creating omegaRD perimeter
sP.mesh       = mesh;
sP.uMesh      = OmegaRD;
sP.epsilon    = mesh.computeMeanCellSize();
sP.filter     = perFilter;
sP.value0     = 1;
P             = PerimeterFunctional(sP);
[POmegaRD,~]    = P.computeFunctionAndGradient(d);