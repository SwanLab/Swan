%% Testing other functions
% Create a Mesh FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

%% Create functions
% AnalyticalFunction

sAF.fHandle = @(x) x(:,2,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);