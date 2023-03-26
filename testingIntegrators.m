%% Testing integrators
% Create a Mesh FEM results
clear; close all;

file = 'test2d_micro';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

%% Create functions
% AnalyticalFunction

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);
p1fun = xFun.project('P1');

% RHS integrator

s.mesh = mesh;
s.type = 'ShapeFunctionFun';
rhs = RHSintegrator.create(s);
integ = rhs.compute(p1fun);