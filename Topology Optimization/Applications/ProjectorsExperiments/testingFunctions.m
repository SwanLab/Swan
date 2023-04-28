%% Testing other functions
% Create a Mesh FEM results
clear; close all;

x =linspace(0,1,2);
y =linspace(0,1,2);

[xv,yv] = meshgrid(x,y);
sM.coord(:,1) = xv(:);
sM.coord(:,2) = yv(:);
sM.connec = delaunay(sM.coord);

sM.coord = [0 0;1 0;0 1];
sM.connec = [1 2 3];

% sM.coord = [0 0;1 0;0 1];
% sM.connec = [2 3 1];

m = Mesh(sM);
mesh = m;

%% Create functions
% AnalyticalFunction

sAF.fHandle = @(x,y) [x(1,:,:),x(2,:,:)];
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

%% Create projectors to P0, P1 and P1D
% % Projector to P0
% pp0.mesh   = mesh;
% pp0.connec = mesh.connec;
% projP0 = Projector_toP0(pp0);
% p0fun = projP0.project(xFun);
% p0fun.plot()
% title('P0')

%% Projector to P1
clc
pp1.mesh   = mesh;
pp1.connec = mesh.connec;
pp1.polynomialOrder = 1;
projP1 = FE_Projector(pp1);
p1fun = projP1.project(xFun);
p1fun.plot()
title('Order1')

% pp1.polynomialOrder = 2;
% projP1 = FE_Projector(pp1);
% p1fun = projP1.project(xFun);
% p1fun.plot()
% title('Order2')
% 
% pp1.polynomialOrder = 3;
% projP1 = FE_Projector(pp1);
% p1fun = projP1.project(xFun);
% p1fun.plot()
% title('Order3')
% 
% pp1.polynomialOrder = 4;
% projP1 = FE_Projector(pp1);
% p1fun = projP1.project(xFun);
% p1fun.plot()
% title('Order4')

% projP12 = Projector_toP1(pp1);
% p1fun2 = projP12.project(xFun);
% p1fun2.plot()
% title('P1 (quad linear) old')

%% Projector to P1 Discontinuous
% pp1d.mesh   = mesh;
% pp1d.connec = mesh.connec;
% projP1D = Projector_toP1Discontinuous(pp1d);
% p1dfun = projP1D.project(xFun);
% p1dfun.plot()
% title('P1 disc (quad quadratic)')

