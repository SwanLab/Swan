clear;
clc;
close all;

filename = 'anisoCantilever';
a.fileName = filename;
gid = FemDataContainer(a);
mesh = gid.mesh;

s.fHandle = @(x) 1-heaviside((x(1,:,:)-1).^2+(x(2,:,:)-0.5).^2-0.3.^2);
s.ndimf   = 1;
s.mesh    = mesh;
fun       = AnalyticalFunction(s);

s.mesh  = mesh;
s.alpha = 4;
s.beta  = 0;
s.theta = 90;
filter  = NonLinearFilterSegment(s);

rhoEps = filter.compute(fun,2);
rhoEps.plot();