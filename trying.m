close all
clear all
mesh = UnitQuadMesh(50,50);
figure()
mesh.plot()

s.fHandle = @(x) [sin(2*pi*x(1,:,:));cos(2*pi*x(1,:,:))]; % f(x) = sin(2*pi*x)
s.mesh    = mesh;
xFun = AnalyticalFunction(s);
plot(project(xFun,'P1'))
