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

% %% Generating a 2D mesh with a hole inclusion
% % Using functions!
% clear; close all
% 
% % Create the data container for the FEM problem
% a.fileName = 'test2d_micro';
% m = FemDataContainer(a);
% 
% % Create the characteristic function (1 inside circle, 0 outside)
% s.mesh    = m.mesh;
% s.fxy     = @(x,y) (x-0.5).^2+(y-0.5).^2-0.3.^2;
% circleFun = CharacteristicFunction(s);