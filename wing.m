clc;clear;close all

addpath(genpath(fileparts(mfilename('fullpath'))))

% Generate coordinates
x1 = linspace(0,30,5);
x2 = linspace(0,2,5);
% Create the grid
[xv,yv] = meshgrid(x1,x2);
% Triangulate the mesh to obtain coordinates and connectivities
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

s.coord = V(:,1:2);
s.connec = F;
mesh = Mesh(s);

figure()
mesh.plot()