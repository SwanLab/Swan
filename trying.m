close all
clear all
mesh = UnitTriangleMesh(50,50);
figure()
mesh.plot()

s.fHandle = @(x) [sin(2*pi*x(1,:,:).*x(2,:,:));cos(2*pi*x(1,:,:).*x(2,:,:))]; % f(x) = sin(2*pi*x)
s.mesh    = mesh;
xFun = AnalyticalFunction(s);
t = xFun.project('P1');
t.plot()
%plot(project(xFun,'P1'))
