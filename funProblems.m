%% Testing funcitons for FEM
clc; clear; close all;

file = 'test2d_triangle';
% file = 'test2d_quad';
% file = 'test3d_hexahedra';
a.fileName = file;
s = FemDataContainer(a);
fem = FunElasticProblem(s);
fem.solve();

%% Not so fast
clear; % close all;

% file = 'test2d_triangle';
file = 'test2d_micro';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

% AnalyticalFunction
sAF.fHandle = @(x) [x(1,:,:).^2; x(2,:,:)];
% sAF.fHandle = @(x) [cos(x(1,:,:).*x(2,:,:)); x(1,:,:).*x(2,:,:)];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

% Quadrature
quad = Quadrature.set(s.mesh.type);
quad.computeQuadrature('LINEAR');


% Projector to P1
pp1.mesh   = mesh;
pp1.connec = mesh.connec;
projP1 = Projector_toP1(pp1);
p1fun = projP1.project(xFun);
p1fun.plot(mesh)

% Gradients
grad1 = p1fun.computeGradient(quad,mesh);
gradientOp = Gradient();
grad2 = gradientOp.compute(p1fun, quad, mesh);

%% Boundary conditions as functions
% AnalyticalFunction
sAF.fHandle = @(x) [x(1,:,:).^2; x(2,:,:)];
% sAF.fHandle = @(x) [cos(x(1,:,:).*x(2,:,:)); x(1,:,:).*x(2,:,:)];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);