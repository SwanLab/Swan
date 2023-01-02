clc; clear; close all;

% %% Export STL from GiD mesh
% file = 'test_micro_holeinclusion';
% % file = 'test2d_triangle';
% a.fileName = file;
% s = FemDataContainer(a);
% fem = FEM.create(s);
% fem.computeChomog();
% % fem.solve();
% fem.print(file)
% 
% P = s.mesh.coord;
% T = s.mesh.connec;
% TR = triangulation(T,P);
% triplot(TR)
% stlwrite(TR, 'whatever.stl') % export as stl

%% Export STL from Unfitted Mesh
clc; clear
% Background

x1 = linspace(-1,1,100);
x2 = linspace(0,1,100);

[xv,yv] = meshgrid(x1,x2);
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');
s.coord  = V(:,1:2);
s.connec = F;
backgroundMesh = Mesh(s);

% Boundary
sB.backgroundMesh = backgroundMesh;
sB.dimension = 1:3;
sB.type = 'FromReactangularBox';
bMc = BoundaryMeshCreator.create(sB);
boundaryMesh  = bMc.create();

% Level set

s.type       = 'circleInclusion';
s.mesh       = backgroundMesh;
s.ndim       = 2;
s.fracRadius = 0.4;
levelSet = LevelSetCreator.create(s);
ls = levelSet.getValue();

% Unfitted mesh
s.boundaryMesh   = boundaryMesh;
s.backgroundMesh = backgroundMesh;
uMesh = UnfittedMesh(s);
uMesh.compute(ls);
uMesh.plot();

% New mesh
% s.coord  = uMesh.innerMesh.mesh.coord;
% s.connec = [uMesh.innerMesh.mesh.connec; uMesh.innerCutMesh.mesh.connec];
% 
% mT = Mesh(s);
% mT.plot

coordInner     = uMesh.innerMesh.mesh.coord;
connecInner    = uMesh.innerMesh.mesh.connec;
coordCutInner  = uMesh.innerCutMesh.mesh.coord;
connecCutInner = uMesh.innerCutMesh.mesh.connec;
ncoord = size(coordInner,1);
connecCutInner = connecCutInner + ncoord;
s.coord = [coordInner;coordCutInner];
s.connec = [connecInner;connecCutInner];
jointMesh = Mesh(s);

% Export to STL
figure()
P = jointMesh.coord;
T = jointMesh.connec;
TR = triangulation(T,P);
triplot(TR)
stlwrite(TR, 'whateverLS.stl')