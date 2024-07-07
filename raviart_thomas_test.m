addpath(genpath(pwd))
close all

% VOLUMES
s.geometryType = "Volume";

% Single tetrahedra
% s.coord = [0,0,0;1,0,0;0,1,0;0,0,2];
% s.connec = [1 2 3 4];

% Unit tetra mesh
mtemp = UnitTetraMesh(5,5,5);
s.coord = mtemp.coord-0.5;
s.connec = mtemp.connec;

m = Mesh.create(s);

sAF.fHandle = @(x) [x(1,:,:); x(2,:,:); x(3,:,:)];
sAF.ndimf   = 3;
sAF.mesh    = m;
xFun = AnalyticalFunction(sAF);


% SURFACES
% s.geometryType = "Surface";

% Single triangle
% s.coord = [0,0;1,0;0,1];
% s.connec = [1 2 3];

% Double triangle
% s.coord = [0,0;1,0;0,1;1,1];
% s.connec = [1 2 3;2 4 3];

% Unit triangle mesh
% mtemp = UnitTriangleMesh(3,3);
% s.coord = mtemp.coord-0.5;
% s.connec = mtemp.connec;
 
% m = Mesh.create(s);
% sAF.fHandle = @(x) [x(2,:,:); -x(1,:,:)];
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