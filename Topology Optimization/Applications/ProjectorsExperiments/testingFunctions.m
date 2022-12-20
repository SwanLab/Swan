%% Testing other functions
% Create a Mesh FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

%% Create functions
% AnalyticalFunction

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

%% Create projectors to P0, P1 and P1D
% Projector to P0
pp0.mesh   = mesh;
pp0.connec = mesh.connec;
projP0 = Projector_toP0(pp0);
resAFtoP0 = projP0.project(xFun);
resAFtoP0.plot(mesh)
title('P0')

% Projector to P1
pp1.mesh   = mesh;
pp1.connec = mesh.connec;
projP1 = Projector_toP1(pp1);
resAFtoP1 = projP1.project(xFun);
resAFtoP1.plot(mesh)
title('P1 (quad linear)')

% Projector to P1 Discontinuous
pp1d.mesh   = mesh;
pp1d.connec = mesh.connec;
projP1D = Projector_toP1Discontinuous(pp1d);
resAFtoP1D = projP1D.project(xFun);
resAFtoP1D.plot(mesh)
title('P1 disc (quad quadratic)')

