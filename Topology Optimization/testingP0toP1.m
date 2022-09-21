%% P0 to P1
% Create sample FEM results
% WIP ton
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

strain = squeeze(fem.variables.strain)';


%% Create FeFunc strain
z.fValues = strain;
strainFun = P0Function(z);

%% Create projector
bb.mesh   = s.mesh;
bb.connec = s.mesh.connec;
projector = Projector_P0toP1(bb);
p1strain = projector.project(strainFun);

%% Use this to create results for testing:
nrows = numel(p1strain.fValues);
xP = reshape(p1strain.fValues,[nrows,1]);
% save('name.mat','xP');