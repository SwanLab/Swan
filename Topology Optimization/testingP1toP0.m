%% P0 to P1
% Create sample FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

uCol   = fem.variables.d_u;
u      = reshape(uCol,[s.mesh.ndim,s.mesh.nnodes])';

%%  Create a FeFunction
z.connec = s.mesh.connec;
z.type   = s.mesh.type;
z.fNodes = u;
uFun = P1Function(z);
uFunD = uFun.computeDiscontinuousField();


%% Create the projector
% cc.mesh   = s.mesh;
% cc.connec = s.mesh.connec;
% cc.nelem  = size(s.mesh.connec,1);
% cc.nnode  = size(s.mesh.connec,2);
% cc.npnod  = size(s.mesh.coord,1);
% projector2 = Projector_P1toP0(cc);
% u_P0 = projector2.project(uFun);

%% Plot using P0Function
% u_P0.plot(s.mesh);