%% Testing gradients from FEM
% Get FEM results

clear; close all;

% file = 'test2d_triangle';
% file = 'test2d_quad';
file = 'Cantileverbeam_Quadrilateral_Bilinear';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();
femU  = reshape(fem.variables.d_u,[s.mesh.ndim,s.mesh.nnodes])';

% P1 Function
z.mesh    = s.mesh.type;
z.fValues = femU;
z.order = 'P1';
uP1 = LagrangianFunction(z);

% Quadrature
quad = Quadrature.set(s.mesh.type);
quad.computeQuadrature('LINEAR');

% Calculate gradient
grad = uP1.computeGradient(quad,s.mesh);
symGrad1 = uP1.evaluateSymmetricGradient(quad.posgp);
symGrad2 = uP1.evaluateSymmetricGradient(quad.posgp);

% Plot
grad.plot(s.mesh);
symGrad1.plot(s.mesh);
symGrad2.plot(s.mesh);

%% Testing gradients via AnalyticalFunction
clear; % close all;

% file = 'test2d_triangle';
file = 'test2d_micro';
% file = 'Cantileverbeam_Quadrilateral_Bilinear';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

% AnalyticalFunction
sAF.fHandle = @(x) x(1,:,:).^2;
sAF.ndimf   = 1;
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
% p1fun.plot(mesh)

symGrad = p1fun.computeSymmetricGradient(quad, mesh);
symGrad.plot(mesh);

% Actual gradient
gradAF = p1fun.computeGradient(quad,mesh);
gradAF.plot(mesh);

%% Testing gradients via AnalyticalFunction
clear; % close all;

file = 'test2d_triangle';
% file = 'test2d_micro';
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

% Actual gradient
gradAF = p1fun.computeGradient(quad,mesh);
gradAF.plot(mesh);

% Symmetric gradient
gradSym1 = p1fun.computeSymmetricGradient(quad, mesh);
gradSym2 = p1fun.computeSymmetricGradient2(quad, mesh);
gradSym1.plot(mesh)
gradSym2.plot(mesh)

