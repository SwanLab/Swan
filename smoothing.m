
clear;
close all;
clc;

load('periodicExample2.mat')

a.fileName = 'SquareForAniTests';
s = FemDataContainer(a);

l.type = 'LevelSet';
l.mesh = s.mesh;
l.creatorSettings.value = ls;
l.initialCase = 'given';
LS = LevelSet(l);

ss.mesh = s.mesh;
ss.fValues = LS.value;
ss.filename = 'seminarPure';
chi = P1Function(ss);
chi.print(ss);


f.mesh = s.mesh;
f.femSettings.scale = 'MACRO';
f.quadratureOrder = 'LINEAR';
f.designVariable = LS;
f.femSettings.LHStype = 'DiffReactRobin';
chiFilter = Filter_PDE_LevelSet(f);
eh = 3;


epsilon = eh*s.mesh.computeMeanCellSize();
chiFilter.updateEpsilon(epsilon);
regularizedDensity = chiFilter.getP1fromP1(LS.value);
ss.fValues = regularizedDensity;
ss.filename = 'seminarExample';
ss.type = 'GiD';
rhoe = P1Function(ss);
rhoe.print(ss);