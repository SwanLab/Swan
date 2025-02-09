addpath(genpath(pwd))
close all

% VOLUMES
s.coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 1 2]+1;
s.connec = [1 2 3 4;5 3 2 4];

% SURFACES
% s.geometryType = "Surface";
% mtemp = UnitTriangleMesh(10,10);
% s.coord = mtemp.coord-0.5;
% s.connec = mtemp.connec;

m = Mesh.create(s);
% sAF.fHandle = @(x) [-x(2,:,:); x(1,:,:)];
sAF.fHandle = @(x) [1-x(2,:,:)-x(3,:,:); x(1,:,:); x(1,:,:)];
sAF.ndimf   = 3;
sAF.mesh    = m;
xFun = AnalyticalFunction(sAF);


% FUNCTIONS
% p1fun = xFun.project('P3');
% rtfun = xFun.project('RT');
nfun  = xFun.project('N');

% p1fun.plot()
% rtfun.plot()
nfun.plot()