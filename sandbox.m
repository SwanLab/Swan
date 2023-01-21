%% This is a sandbox file!
% Feel free to test anything here :)
clc; clear; close all;

% file = 'test2d_triangle';
% file = 'test2d_quad';
file = 'test3d_hexahedra';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;
fem = FEM.create(s);
fem.solve();

%% Create functions
% AnalyticalFunction

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
% sAF.fHandle = @(x) [x(1,:,:).^2; x(2,:,:)];
% sAF.ndimf   = 2;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

%% Create projectors to P0, P1 and P1D
% testingFunctions.m
% testingGradients.m
% Projector to P1
pp1.mesh   = mesh;
pp1.connec = mesh.connec;
projP1 = Projector_toP1(pp1);
p1fun = projP1.project(xFun);
% p1fun.plot(mesh)
% title('P1 (quad linear)')

%% Function printing
aa.mesh = mesh;
aa.filename = 'hellothere';
p1fun.print(aa)