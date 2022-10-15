%% Testing gradients from FEM
% Get FEM results
clear; close all;

file = 'test2d_quad';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();
femU  = reshape(fem.variables.d_u,[s.mesh.ndim,s.mesh.nnodes])';

% P1 Function
z.connec  = s.mesh.connec;
z.type    = s.mesh.type;
z.fValues = femU;
uP1 = P1Function(z);

% Quadrature
quad = Quadrature.set(s.mesh.type);
quad.computeQuadrature('LINEAR');

% Calculate gradient
fg = uP1.computeSymmetricGradient(quad, s.mesh);

% Project
pp1.mesh   = s.mesh;
pp1.connec = s.mesh.connec;
projP1 = Projector_toP1(pp1);
p1strain = projP1.project(fg);
p1strain.plot(s.mesh);

%% Testing gradients via AnalyticalFunction
clear; % close all;

% file = 'test2d_triangle';
file = 'test2d_micro';
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

fgauss = p1fun.computeSymmetricGradient(quad, mesh);
p1fg = projP1.project(fgauss);
p1fg.plot(mesh);

% Actual gradient
gradAF = p1fun.computeGradient(quad,mesh);
gradAF.plot(mesh);

%% Testing gradients via AnalyticalFunction
clear; % close all;

% file = 'test2d_triangle';
file = 'test2d_micro';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;

% AnalyticalFunction
sAF.fHandle = @(x) [x(1,:,:).^2; x(2,:,:)];
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
% p1fun.plot(mesh)

fgauss = p1fun.computeSymmetricGradient(quad, mesh);
p1fg = projP1.project(fgauss);
% p1fg.plot(mesh);

% Actual gradient
gradAF = p1fun.computeGradient(quad,mesh);
gradAF.plot(mesh);
