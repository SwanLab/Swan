%% Testing projection to P0
% Get FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

femStrain = squeeze(fem.variables.strain)';
femU      = reshape(fem.variables.d_u,[s.mesh.ndim,s.mesh.nnodes])';

%% Create functions
% P0 Function
% z.fValues = femStrain;
z.fValues = femStrain(:,1:2);
% z.fValues = femStrain(:,2);
z.mesh    = s.mesh;
z.order   = 'P0';
sigP0 = LagrangianFunction(z);

% P1 Function
z.mesh    = s.mesh;
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
quad = Quadrature.set(s.mesh.type);
quad.computeQuadrature('LINEAR');
x.fValues    = permute(fem.variables.strain,[2 1 3]);
x.connec     = s.mesh.connec;
x.type       = s.mesh.type;
x.quadrature = quad;
sigFg = FGaussDiscontinuousFunction(x);

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
resFgtoP0 = projP0.project(sigFg);

%% Projection to P1
% Create the projector
pp1.mesh   = s.mesh;
pp1.connec = s.mesh.connec;
projP1 = Projector_toP1(pp1);

% P0 to P1
resP0toP1 = projP1.project(sigP0);

% P1 Discontinuous to P1
resP1DtoP1 = projP1.project(uP1D);

% FGauss function to P1
resFgtoP1 = projP1.project(sigFg);

%% Projection to P1 Discontinuous
% Create the projector
pp1d.mesh   = s.mesh;
pp1d.connec = s.mesh.connec;
projP1D = Projector_toP1Discontinuous(pp1d);

% P0 to P1 Discontinuous
resP0toP1D = projP1D.project(sigP0);

% P1 to P1 Discontinuous
resP1toP1D = projP1D.project(uP1);
resP1toP1D.plot(s.mesh)

% FGauss function to P1 Discontinuous
resFgtoP1D = projP1D.project(sigFg);