% Checking double contraction
clc; clear; close all;

mesh = UnitQuadMesh(20,20);
gPar.type          = 'Circle';
gPar.radius        = 0.25;
gPar.xCoorCenter   = 0.5;
gPar.yCoorCenter   = 0.5;
g                  = GeometricalFunction(gPar);
phiFun             = g.computeLevelSetFunction(mesh);
lsCircle           = phiFun.fValues;
lsCircleInclusion  = -lsCircle;
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(lsCircleInclusion);


%% Create Inner Mesh
% ONLY using MATLAB (GiD improves conditioning)
IM = uMesh.createInnerMesh();
IM.plot