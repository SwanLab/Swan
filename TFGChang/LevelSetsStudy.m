clear;
close all;

Naca.M     = 0.02;
Naca.p     = 0.4;
Naca.t     = 0.12;
Naca.chord = 1;
Naca.AoA   = 0;

length  = 8;
height  = 4;
chord   = 1;
nx      = 600;
ny      = 300;
refMesh = TriangleMesh(length,height,nx,ny);


% Example LS1
s.type = 'LevelSet1';
s.xLE  = (length - chord)/2;
s.yLE  = height/2;

s.chord = Naca.chord;
s.p     = Naca.p;
s.m     = Naca.M;
s.t     = Naca.t;
s.AoA   = Naca.AoA;

g  = GeometricalFunction(s);
levelSet1 = g.computeLevelSetFunction(refMesh);
levelSet1.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet1.fValues);
m1 = uMesh.createInnerMesh();
figure
m1.plot();



% Example LS2
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
levelSet2.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(levelSet2.fValues);
m2 = uMesh.createInnerMesh();
figure
m2.plot();


ls = exp(20*levelSet1.*levelSet2) - 1;
ls = ls.project('P1');
lsVal = ls.fValues;
lsVal(abs(lsVal)<=1e-4) = -1e-4; % DELETE OSCILLATIONS WITH L1 FILTER?
ls.setFValues(lsVal);
ls.plot();

s.backgroundMesh = refMesh;
s.boundaryMesh   = refMesh.createBoundaryMesh();
uMesh        = UnfittedMesh(s);
uMesh.compute(ls.fValues);
m3 = uMesh.createInnerMesh();
figure
m3.plot();