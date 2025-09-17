close all
clear all
mesh = UnitQuadMesh(50,50);
figure()
mesh.plot()

sAF.fHandle = @(x) sin(2*pi*x(1,:,:)); % f(x) = sin(2*pi*x)
sAF.ndimf   = 1; % number of dimensions
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);
plot(project(xFun,'P1'))
