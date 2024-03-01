close all
clear
clc

s.coord  = [-1 -1 -1;
            +1 -1 -1;
            +1 +1 -1;
            -1 +1 -1;
            -1 -1 +1;
            +1 -1 +1;
            +1 +1 +1;
            -1 +1 +1;
            -1 -1 +3;
            +1 -1 +3;
            +1 +1 +3;
            -1 +1 +3];

s.connec = [1 2 3 4 5 6 7 8;
    5 6 7 8 9 10 11 12];

% s.coord  = [-1 -1 -1;
%             +1 -1 -1;
%             -1 +1 -1;
%             -1 -1 +1;
%             +1 +1 +1];
% 
% s.connec = [1 2 3 4;
%             2 3 4 5];

mesh = Mesh.create(s);

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1fun = xFun.project('P1');
p2fun = xFun.project('P2');
p3fun = xFun.project('P3');
