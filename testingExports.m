%% Testing mesh exporting
clc; clear; close all;

% Create background mesh
x1 = linspace(-1,1,20);
x2 = linspace(0,1,20);
[xv,yv] = meshgrid(x1,x2);
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
sBg.coord  = V(:,1:2);
sBg.connec = F;
bgMesh = Mesh(sBg);

% Create boundary mesh
sBd.backgroundMesh = bgMesh;
sBd.dimension      = 1:3;
sBd.type = 'FromReactangularBox';
bMc = BoundaryMeshCreator.create(sBd);
bdMesh  = bMc.create();

% Create Level Set
sLS.type       = 'circleInclusion';
sLS.mesh       = bgMesh;
sLS.ndim       = 2;
sLS.fracRadius = 0.4;
sLS.coord      = bgMesh.coord;
ls = LevelSetCreator.create(sLS);
levelSet = ls.getValue();

% Create Unfitted Mesh
sUm.boundaryMesh   = bdMesh;
sUm.backgroundMesh = bgMesh;
uMesh = UnfittedMesh(sUm);
uMesh.compute(levelSet);

%% Create Inner Mesh
% ONLY using MATLAB (GiD improves conditioning)
IM = uMesh.createInnerMesh();

%% Create Inner Mesh And Improve Conditioning
% ONLY using GiD (MATLAB does not improve conditioning)
IMcond = uMesh.createInnerMeshGoodConditioning();

%% Extrude Mesh
height = 0.16;
EM = IMcond.provideExtrudedMesh(height);

%% Export STL
EM.exportSTL();