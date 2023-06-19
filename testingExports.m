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

%% 
% IM.improveConditioning(); % only for 2d, using 3d


%% Create Inner Mesh And Improve Conditioning
% ONLY using GiD (MATLAB does not improve conditioning)
sIMg.filename     = 'notForLong';  % Why?
sIMg.meshFileName = 'notForLong2'; % Why?
sIMg.swanPath     = '/home/ton/Github/Swan/';
sIMg.gidPath      = '/home/ton/GiDx64/gid-16.1.2d/';
% IMcond = uMesh.createInnerMeshGoodConditioning(sIMg);
IMcond2 = uMesh.createInnerMeshGoodConditioning(sIMg);

%% Extrude Mesh (improves conditioning)
sEM.filename     = 'notForLong';  % Why?
sEM.meshFileName = 'notForLong2'; % Why?
sEM.swanPath     = '/home/ton/Github/Swan/';
sEM.gidPath      = '/home/ton/GiDx64/gid-16.1.2d/';
EM = uMesh.provideExtrudedMesh(sEM);

%% Export STL
sSTL.filename     = 'notForLong';  % Why?
sSTL.meshFileName = 'notForLong2'; % Why?

EM.exportSTL(sSTL); % 2 options: MATLAB/Gid
% Matlab 1 -> export stl using tetrahedra -> searhco nline
% Matlab 2 -> use boundary mesh and boundary cut mesh...

%Mescape utilities SwapNormals Surfaces MakeGroupCoherent 1:END escape Yes escape


% In short: research STL