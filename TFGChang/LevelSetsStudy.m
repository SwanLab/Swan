clear;
%close all;

% Naca Info
Naca.M     = 0.02;
Naca.p     = 0.4;
Naca.t     = 0.12;
Naca.chord = 1;
Naca.AoA   = 0;

% Mesh Params
length  = 8;
height  = 4;
chord   = 1;
nx      = 600;
ny      = nx/2;
refMesh = TriangleMesh(length,height,nx,ny);

% %% Example Test
% % clear;
% % close all;
% s.type = 'LevelSetTest';
% s.xLE  = (length - chord)/2;
% s.yLE  = height/2;
% 
% s.chord = Naca.chord;
% s.p     = Naca.p;
% s.m     = Naca.M;
% s.t     = Naca.t;
% s.AoA   = Naca.AoA;
% 
% g  = GeometricalFunction(s);
% levelSetTest = g.computeLevelSetFunction(refMesh);
% levelSetTest.plot();
% 
% s.backgroundMesh = refMesh;
% s.boundaryMesh   = refMesh.createBoundaryMesh();
% uMesh        = UnfittedMesh(s);
% uMesh.compute(levelSetTest.fValues);
% mT = uMesh.createInnerMesh();
% figure
% mT.plot();


%% Example LS1
s.type = 'LevelSet1';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSetTest = g.computeLevelSetFunction(refMesh);
%levelSetTest.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSetTest.fValues);
m1 = uMesh.createInnerMesh();
figure
m1.plot();


%% Example LS2
s.type = 'LevelSet2';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet2 = g.computeLevelSetFunction(refMesh);
%levelSet2.plot();


s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet2.fValues);
m2 = uMesh.createInnerMesh();
figure
m2.plot();

m = collapseMeshes(m1,m2);

points = m.coord;
r = inf;
T      = alphaShape(points,r);
DT     = alphaTriangulation(T);

s.connec = DT;
s.coord  = points;
m        = Mesh.create(s);

s.type = 'LevelSet1';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet1 = g.computeLevelSetFunction(m);

s.type = 'LevelSet2';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet2 = g.computeLevelSetFunction(m);

q        = Quadrature.create(m, 0);
xV       = q.posgp;

lsElem1 = squeeze(levelSet1.evaluate(xV));
lsElem2 = squeeze(levelSet2.evaluate(xV));
s.connec = m.connec(lsElem1.*lsElem2<=0,:);

s.coord        = m.coord;
m2             = Mesh.create(s);

ls = exp(20*levelSet1.*levelSet2) - 1;
ls = ls.project('P1');
lsVal = ls.fValues;
lsVal(abs(lsVal)<=1e-4) = -1e-4; % DELETE OSCILLATIONS WITH L1 FILTER?
ls.setFValues(lsVal);
ls.plot();
m2.plot();

%% Example LS3
s.type = 'LevelSet3';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet3 = g.computeLevelSetFunction(refMesh);
levelSet3.plot();


s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(ls.fValues);
m3 = uMesh.createInnerMesh();
figure
m3.plot();



%% Example LS4
s.type = 'LevelSet4';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet4 = g.computeLevelSetFunction(refMesh);
levelSet4.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet4.fValues);
m4 = uMesh.createInnerMesh();
% figure
m4.plot();


%% Example LS4
s.type = 'NacaHole';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet4 = g.computeLevelSetFunction(refMesh);
levelSet4.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet4.fValues);
m4 = uMesh.createInnerMesh();
 figure
m4.plot();

%% Example yc
s.type = 'LevelSet5';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet5 = g.computeLevelSetFunction(refMesh);
levelSet5.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet5.fValues);
m5 = uMesh.createInnerMesh();
 figure
m5.plot();

%% Example yt
s.type = 'LevelSet6';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet6 = g.computeLevelSetFunction(refMesh);
 levelSet6.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet6.fValues);
m6 = uMesh.createInnerMesh();
 figure
m6.plot();

%% Example yu
s.type = 'LevelSet7';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet7 = g.computeLevelSetFunction(refMesh);
levelSet7.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet7.fValues);
m7 = uMesh.createInnerMesh();
figure
m7.plot();

