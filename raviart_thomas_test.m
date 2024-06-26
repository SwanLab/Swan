addpath(genpath(pwd))
close all

s.geometryType = "Surface";

s.coord = [0,0;1,0;0,1;1,1];
s.connec = [1 2 3; 4 2 3];

s.coord = [0,0;1,0;0,1];
s.connec = [1 2 3];

mtemp = UnitTriangleMesh(5,5);
s.coord = mtemp.coord-0.5;
s.connec = mtemp.connec;
m = Mesh.create(s);

rt = RaviartThomasFunction.create(m,1,1);

sAF.fHandle = @(x) [x(1,:,:) ; x(2,:,:)];
sAF.ndimf   = 2;
sAF.mesh    = m;
xFun = AnalyticalFunction(sAF);

p1fun = xFun.project('P1');
rtfun = xFun.project('RT');
p1fun.plot()
rtfun.plot()