close all
clear
clc

s.coord  = [0;1;2];
s.connec = [1 2;2 3];

mesh = Mesh.create(s);

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1fun = xFun.project('P1');
p2fun = xFun.project('P2');
p3fun = xFun.project('P3');