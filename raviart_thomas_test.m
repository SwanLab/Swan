addpath(genpath(pwd))
close all

% VOLUMES
s.geometryType = "Volume";
mtemp = UnitTetraMesh(2,2,2);
s.coord = mtemp.coord-0.5;
s.connec = mtemp.connec;

m = Mesh.create(s);

sAF.fHandle = @(x) [x(1,:,:); x(2,:,:); x(3,:,:)];
sAF.ndimf   = 3;
sAF.mesh    = m;
xFun = AnalyticalFunction(sAF);


% SURFACES
% s.geometryType = "Surface";
% mtemp = UnitTriangleMesh(2,2);
% s.coord = mtemp.coord-0.5;
% s.connec = mtemp.connec;
% 
% m = Mesh.create(s);
% sAF.fHandle = @(x) [x(1,:,:); x(2,:,:)];
% sAF.ndimf   = 2;
% sAF.mesh    = m;
% xFun = AnalyticalFunction(sAF);


% FUNCTIONS
% p1fun = xFun.project('P1');
rtfun = xFun.project('RT');
% nfun  = xFun.project('N');

% p1fun.plot()
rtfun.plot()
% nfun.plot()