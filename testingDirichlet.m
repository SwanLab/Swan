%% Testing subfunctions
% Goal: generate subfunctions of functions (eg. left boundary)
clear; close all;

% Create P1 Function
sAF.fHandle = @(x) [x(1,:,:).^2];
sAF.ndimf   = 1;
sAF.mesh    = UnitQuadMesh(7,7);
xFun = AnalyticalFunction(sAF);

p1 = xFun.project('P1');