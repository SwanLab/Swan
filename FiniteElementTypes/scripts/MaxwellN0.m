%% Maxwell equations with first order Nedelec elements
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

% sM.coord = [0 1;0 0;1 0];
% sM.connec = [1 2 3];

% sM.coord = [1 1;0 1;1 0];
% sM.connec = [1 2 3];

mesh = Mesh(sM);

%% Create functions
% AnalyticalFunction

sAF.fHandle = @(x,y) [x(1,:,:),x(2,:,:)];
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

%% Projector to P1
clc
pp1.mesh   = mesh;
pp1.connec = mesh.connec;
pp1.polynomialOrder = 1;
pp1.feParams.type = "Nedelec";
pp1.feParams.order = 1;
pp1.feParams.dim = 2;
projP1 = FE_Projector(pp1);
p1fun = projP1.project(xFun);
p1fun.plot()
title('Order1')