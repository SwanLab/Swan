close all;
clear;
clc

mesh = UnitQuadMesh(1,1);

sAF.fHandle = @(x) sin(2*pi*x(1,:,:));
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1fun = xFun.project('P2');
% p2fun = xFun.project('P2');
% p3fun = xFun.project('P3');