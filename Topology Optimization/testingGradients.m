%% Testing gradients
% Get FEM results
clear; close all;

file = 'test2d_quad';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

femU      = reshape(fem.variables.d_u,[s.mesh.ndim,s.mesh.nnodes])';

%% Create functions

% P1 Function
z.connec  = s.mesh.connec;
z.type    = s.mesh.type;
z.fValues = femU;
uP1 = P1Function(z);

% Quadrature
quad = Quadrature.set(s.mesh.type);
quad.computeQuadrature('LINEAR');

% Calculate gradient
fg = uP1.computeGradientStrain(quad, s.mesh);

% Project
pp1.mesh   = s.mesh;
pp1.connec = s.mesh.connec;
projP1 = Projector_toP1(pp1);
p1strain = projP1.project(fg);
