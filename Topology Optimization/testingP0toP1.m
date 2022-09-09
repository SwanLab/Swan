%% Create sample FEM results
% WIP ton
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

strain = squeeze(fem.variables.strain)';


%% P0 to P1
% Create FeFunc strain
z.mesh    = s.mesh;
z.fValues = strain(:,1,:);
strainFeFun = P0Function(z);


bb.mesh   = s.mesh;
bb.connec = s.mesh.connec;
bb.nelem  = size(s.mesh.connec,1);
bb.nnode  = size(s.mesh.connec,2);
bb.npnod  = size(s.mesh.coord,1);
projector = Projector_P0toP1(bb);


strainCol = reshape(strain, [bb.nelem*3, 1]);
strainP1 = projector.project(strainFeFun);
