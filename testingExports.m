%% Testing mesh exporting
clc; clear; close all;

mesh = UnitQuadMesh(7,7);
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

% %% Create Inner Mesh And Improve Conditioning
% % ONLY using GiD (MATLAB does not improve conditioning)
% IMcond = uMesh.createInnerMeshGoodConditioning();
% 
%% Extrude Mesh
height = 0.16;
% EM = IMcond.provideExtrudedMesh(height);
EM = IM.provideExtrudedMesh(height);

%% Export STL
EM.exportSTL();