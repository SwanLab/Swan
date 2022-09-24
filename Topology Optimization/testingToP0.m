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
z.fValues = femStrain;
% z.fValues = femStrain(:,1:2);
z.connec  = s.mesh.connec;
z.type    = s.mesh.type;
sigP0 = P0Function(z);

% P1 Function
z.connec  = s.mesh.connec;
z.type    = s.mesh.type;
z.fValues = femU;
uP1 = P1Function(z);

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
uFg = FGaussDiscontinuousFunction(x);

%% Projection to P0
% Create the projector
pp.mesh   = s.mesh;
pp.connec = s.mesh.connec;
projP0 = Projector_P1toP0(pp);

% P1 to P0
resP1 = projP0.project(uP1);

% P1 Discontinuous to P0
resP1D = projP0.project(uP1D);

% P1 Discontinuous to P0
resFg = projP0.project(uFg);