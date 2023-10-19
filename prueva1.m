clc;clear; close all

addpath(genpath(fileparts(mfilename('fullpath'))))

x = linspace(0,1,15);
x1 = 0.5 * x;
x2 = 0.5 * x + 2;

plot(x,x1,x,x2)

[xv,yv] = meshgrid(x1,x2);
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

s.coord = V(:,1:2);
s.connec = F;
mesh = Mesh(s);

figure()
mesh.plot()