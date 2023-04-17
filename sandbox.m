%% This is a sandbox file!
% Feel free to test anything here :)
% clc; clear; close all;

% Load mesh
file = 'test2d_micro';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;
clear s;

% LHS integrator
trial = P1Function.create(mesh, 1);
test  = P0Function.create(mesh, 1);

s.type = 'MassTestTrial';
s.mesh = mesh;
s.test = test;
s.trial = trial;
lhs = LHSintegrator.create(s);
LHS = lhs.compute();
