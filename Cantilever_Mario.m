%clc;clear;close all

addpath(genpath(fileparts(mfilename('fullpath'))))

% Generate coordinates
x1 = linspace(0,2,20);
x2 = linspace(1,2,20);
% Create the grid
[xv,yv] = meshgrid(x1,x2);
% Triangulate the mesh to obtain coordinates and connectivities
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

s.coord = V(:,1:2);
s.connec = F;
mesh = Mesh(s);

figure(1)
mesh.plot()
