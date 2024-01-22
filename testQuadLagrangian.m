close all;
clear;
clc

mesh = UnitQuadMesh(2,2);

sAF.fHandle = @(x) x(1,:,:).*x(1,:,:)+x(2,:,:).*x(2,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1fun = xFun.project('P1');
p2fun = xFun.project('P2');
p3fun = xFun.project('P3');