close all;
clear;
clc

s.coord  = [0 0 0;
            0 0 1;
            0 1 0;
            1 0 0;];
s.connec = [1 2 3 4];
mesh = Mesh(s);

sAF.fHandle = @(x) x(1,:,:)+x(2,:,:)+x(3,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

% p1fun = xFun.project('P1');
% p2fun = xFun.project('P2');
% p3fun = xFun.project('P3');
