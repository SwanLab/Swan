%% Testing Operators with functional programming
% Create the Lagrangian Function
clc; clear;

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = UnitTriangleMesh(5,5);
xFun = AnalyticalFunction(sAF);
p1 = xFun.project('P1');

%