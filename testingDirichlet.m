%% Testing subfunctions
% Goal: generate subfunctions of functions (eg. left boundary)
clear; close all;

% Create P1 Function
mesh = UnitQuadMesh(10,10);
sAF.fHandle = @(x) x(2,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

f_omega = xFun.project('P1');

% isMiddle = @(coor) (abs(coor(:,2) >= 0.4 & abs(coor(:,2)) <= 0.6) );
isMiddleRight = @(coor) (abs(coor(:,1)) == 1 & abs(coor(:,2) >= 0.4 & abs(coor(:,2)) <= 0.6) );

% Boundary submesh
[mesh_bound, l2g] = mesh.getBoundarySubmesh(isMiddleRight);
% m.plot

p1_bound = f_omega.restrictTo(isMiddleRight);
dLambda = P0Function.create(mesh_bound, 2);

% LHS
s.test  = dLambda;
s.trial = p1_bound;
s.mesh  = mesh_bound;
s.type  = 'MassMatrix';
lhs = LHSintegrator.create(s);
LHS = lhs.compute();
