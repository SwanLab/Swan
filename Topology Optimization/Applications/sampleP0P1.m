%% P0-P1 Mass matrix example
clc; clear; close all;

% Load mesh
file = 'test2d_micro';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;
clear s;

% Create functions

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1trial = xFun.project('P1');
% p1trial = P1Function.create(mesh, 1);
p0test  = P0Function.create(mesh, 1);

% LHS integrator

s.type = 'MassMatrix';
s.mesh = mesh;
s.test = p0test;
s.trial = p1trial;
lhs = LHSintegrator.create(s);
LHS = lhs.compute();

% Mass P0

s.type = 'MassMatrix';
s.mesh = mesh;
s.test = p0test;
s.trial = p0test;
mp0 = LHSintegrator.create(s);
MP0 = mp0.compute();

% Result

gj = MP0\(LHS*p1trial.fValues);

% Plot
z.fValues = gj;
z.mesh = mesh;
p0_result = P0Function(z);
p0_result.plot();