%% Testing subfunctions
% Goal: generate subfunctions of functions (eg. left boundary)
clear; close all;

% Create P1 Function
mesh = UnitQuadMesh(10,10);
sAF.fHandle = @(x) x(2,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1 = xFun.project('P1');

% Boundary submesh
isMiddle = @(coor) (abs(coor(:,2) >= 0.4 & abs(coor(:,2)) <= 0.6) );
isMiddleRight = @(coor) (abs(coor(:,1)) == 1 & abs(coor(:,2) >= 0.4 & abs(coor(:,2)) <= 0.6) );
% [m, l2g] = mesh.getBoundarySubmesh(isMiddleRight);
% m.plot
p1sub = p1.evaluateBoundarySubdomain(isMiddleRight);

%% Let's do it in P0
clear; close all;

mesh = UnitQuadMesh(10,10);
sAF.fHandle = @(x) x(2,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p0 = xFun.project('P0');
isMiddleRight = @(coor) (abs(coor(:,1)) == 1 & abs(coor(:,2) >= 0.4 & abs(coor(:,2)) <= 0.6) );

% p0sub = p0.evaluateBoundarySubdomain(isMiddleRight);s

[mesh_sub, l2g] = mesh.getBoundarySubmesh(isMiddleRight);