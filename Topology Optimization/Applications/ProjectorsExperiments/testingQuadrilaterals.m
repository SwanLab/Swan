%% Testing quadrilaterals
% Get FEM results
clear; close all;

file = 'test2d_quad';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

femU      = reshape(fem.variables.d_u,[s.mesh.ndim,s.mesh.nnodes])';

%% Create functions
% P0 Function

% It makes no sense -- strain is a FGaussFunction!

% P1 Function
z.connec  = s.mesh.connec;
z.type    = s.mesh.type;
z.fValues = femU;
z.order   = 'P1';
uP1 = LagrangianFunction(z);

% P1 Discontinuous Function
cc.mesh   = s.mesh;
cc.connec = s.mesh.connec;
projectorDisc = ProjectorToP1discont(cc);
ss.origin = 'P1';
ss.x      = uP1;
uP1D = projectorDisc.project(ss);

% FGauss function
quadrature = Quadrature.set(s.mesh.type);
quadrature.computeQuadrature('LINEAR');
x.fValues    = permute(fem.variables.strain,[2 1 3]);
x.connec     = s.mesh.connec;
x.type       = s.mesh.type;
x.quadrature = quadrature;
uFg = FGaussDiscontinuousFunction(x);

%% Projection to P0
% Create the projector
pp0.mesh   = s.mesh;
pp0.connec = s.mesh.connec;
projP0 = Projector_toP0(pp0);

% P1 to P0
resP1toP0 = projP0.project(uP1);

% P1 Discontinuous to P0
resP1DtoP0 = projP0.project(uP1D);

% FGauss function to P0
resFgtoP0 = projP0.project(uFg);

%% Projection to P1
% Create the projector
pp1.mesh   = s.mesh;
pp1.connec = s.mesh.connec;
projP1 = Projector_toP1(pp1);

% % P0 to P1
% resP0toP1 = projP1.project(sigP0);

% P1 Discontinuous to P1
resP1DtoP1 = projP1.project(uP1D);

% FGauss function to P1
resFgtoP1 = projP1.project(uFg);

%% Projection to P1 Discontinuous
% Create the projector
pp1d.mesh   = s.mesh;
pp1d.connec = s.mesh.connec;
projP1D = Projector_toP1Discontinuous(pp1d);

% % P0 to P1 Discontinuous
% resP0toP1D = projP1D.project(sigP0);

% P1 to P1 Discontinuous
resP1toP1D = projP1D.project(uP1);
resP1toP1D.plot(s.mesh)

% FGauss function to P1 Discontinuous
resFgtoP1D = projP1D.project(uFg);

%% UTILS | CREATE A LINEAR QUADRATURE
% quadrature = Quadrature.set(s.mesh.type);
% quadrature.computeQuadrature('LINEAR');
