clear;
close all;

% Naca Info
Naca.M     = 0.09;
Naca.p     = 0.8;
Naca.t     = 0.05;
Naca.chord = 1;
Naca.AoA   = 0;

% Mesh Params
length  = 8;
height  = 4;
chord   = 1;
nx      = 600;
ny      = nx * height/length /0.8;
refMesh = TriangleMesh(length,height,nx,ny);

%% Test Example Using Max of 4 Level Set Functions

s.type = 'LevelSetTest';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSetTest = g.computeLevelSetFunction(refMesh);
levelSetTest.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSetTest.fValues);
mT = uMesh.createInnerMesh();
figure
mT.plot();

%% Test Example Using Collapse Mesh of yu and yl
% yl
s.type = 'LevelSetYU';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet1 = g.computeLevelSetFunction(refMesh);
%levelSet1.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet1.fValues);
m1 = uMesh.createInnerMesh();
% figure
% m1.plot();

% yu
s.type = 'LevelSetYL';

g  = GeometricalFunction(s);
levelSet2 = g.computeLevelSetFunction(refMesh);
%levelSet2.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet2.fValues);
m2 = uMesh.createInnerMesh();
% figure
% m2.plot();

m = collapseMeshes(m2,m1);
figure;
m.plot();

%% Attempt to remove the peak at the leading edge Using GoodMeshConditioning

points = m.coord;
r = inf;
T      = alphaShape(points,r);
DT     = alphaTriangulation(T);

s.connec = DT;
s.coord  = points;
m        = Mesh.create(s);

q        = Quadrature.create(m, 0);
xV       = q.posgp;

lsElem1 = squeeze(levelSet1.evaluate(xV));
lsElem2 = squeeze(levelSet2.evaluate(xV));
s.connec = m.connec(lsElem1.*lsElem2<=0,:);

s.coord        = m.coord;
m2             = Mesh.create(s);

figure;
m.plot();
