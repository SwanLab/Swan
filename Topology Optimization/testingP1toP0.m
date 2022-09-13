%% Create sample FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

uCol   = fem.variables.d_u;
u      = reshape(uCol,[s.mesh.ndim,s.mesh.nnodes])';

%% P1 to P0 v2

% Create a FeFunction
z.connec = s.mesh.connec;
z.type   = s.mesh.type;
z.fNodes = u;
uFeFun = P1Function(z);

% Create the projector
cc.mesh   = s.mesh;
cc.connec = s.mesh.connec;
cc.nelem  = size(s.mesh.connec,1);
cc.nnode  = size(s.mesh.connec,2);
cc.npnod  = size(s.mesh.coord,1);
projector2 = Projector_P1toP0(cc);
u_P0 = projector2.project(uFeFun);

% Plot using P0Function
u_P0.plot(s.mesh);