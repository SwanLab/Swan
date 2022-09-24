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
z.connec  = s.mesh.connec;
z.type    = s.mesh.type;
z.fValues = u;
uFun = P1Function(z);

%% Create the projector
% To P0 
cc.mesh   = s.mesh;
cc.connec = s.mesh.connec;
projector2 = Projector_toP0(cc);
u_P0 = projector2.project(uFun);

% To P1 Discontinuous
projectorDisc = ProjectorToP1discont(cc);
ss.origin = 'P1';
ss.x      = uFun;
uFundisc = projectorDisc.project(ss);

%% Plot using P0Function
% u_P0.plot(s.mesh);

%% Use this to create results for testing:
% nrows = numel(u_P0.fValues);
% xP = reshape(u_P0.fValues,[nrows,1]);
% save('name.mat','xP');
