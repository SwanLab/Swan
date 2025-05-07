clear;
close all;
clc;

% Domain and design variable
load('NumericalExperiments/Gripping/DesVarDensityMinNSACompliance.mat','d');
mesh = d.fun.mesh;

epsOverh = 8;

% Filter PDE for the perimeter
h           = mesh.computeMeanCellSize();
sfi.trial   = LagrangianFunction.create(mesh,1,'P1');
sfi.mesh    = mesh;
sfi.LHStype = 'StiffnessMass';
perFilter   = FilterPDE(sfi);

% Creating local domain hinge A
sG.type            = 'Circle';
sG.xCoorCenter     = 0.28;
sG.yCoorCenter     = 0.5;
sG.radius          = 0.035;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaA             = UnfittedMesh(sUm);
OmegaA.compute(lsFun.fValues);

% Creating local domain hinge B
sG.type            = 'Circle';
sG.xCoorCenter     = 0.475;
sG.yCoorCenter     = 0.745;
sG.radius          = 0.035;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaB             = UnfittedMesh(sUm);
OmegaB.compute(lsFun.fValues);

% Creating local domain hinge C
sG.type            = 'Circle';
sG.xCoorCenter     = 0.475;
sG.yCoorCenter     = 0.255;
sG.radius          = 0.035;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaC            = UnfittedMesh(sUm);
OmegaC.compute(lsFun.fValues);

% Creating local hinge D
sG.type            = 'Circle';
sG.xCoorCenter     = 0.75;
sG.yCoorCenter     = 0.5;
sG.radius          = 0.035;
g                  = GeometricalFunction(sG);
lsFun              = g.computeLevelSetFunction(mesh);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh();
OmegaD            = UnfittedMesh(sUm);
OmegaD.compute(lsFun.fValues);

% Creating omegaA perimeter
sP.mesh      = mesh;
sP.uMesh     = OmegaA;
sP.epsilon   = epsOverh*mesh.computeMeanCellSize();
sP.filter    = perFilter;
sP.value0    = 1;
P            = PerimeterFunctional(sP);
[POmegaA,~] = P.computeFunctionAndGradient(d);

% Creating omegaB perimeter
sP.mesh      = mesh;
sP.uMesh     = OmegaB;
sP.epsilon   = epsOverh*mesh.computeMeanCellSize();
sP.filter    = perFilter;
sP.value0    = 1;
P            = PerimeterFunctional(sP);
[POmegaB,~] = P.computeFunctionAndGradient(d);

% Creating omegaC perimeter
sP.mesh      = mesh;
sP.uMesh     = OmegaC;
sP.epsilon   = epsOverh*mesh.computeMeanCellSize();
sP.filter    = perFilter;
sP.value0    = 1;
P            = PerimeterFunctional(sP);
[POmegaC,~] = P.computeFunctionAndGradient(d);

% Creating omegaD perimeter
sP.mesh      = mesh;
sP.uMesh     = OmegaD;
sP.epsilon   = epsOverh*mesh.computeMeanCellSize();
sP.filter    = perFilter;
sP.value0    = 1;
P            = PerimeterFunctional(sP);
[POmegaD,~] = P.computeFunctionAndGradient(d);