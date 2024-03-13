%% Swan: first steps
% 15th March, 2024
clc; clear; close all

%% Functions
% First we need to generate the mesh
mesh = UnitTriangleMesh(15,15);

% Now we can create the analytical function

sAF.fHandle = @(x) sin(6*x(1,:,:)).*cos(3*x(2,:,:));
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);


%% Projection
% Related classes: Projector_toLagrangian, LHSintegrator, RHSintegrator

p0 = xFun.project('P0');
p1 = xFun.project('P1');
p1.plot

%% Elasticity
% See: ElasticProblem, Tutorial02FEMElasticity