%% Create sample FEM results
clc; clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

strain = fem.variables.strain;
u = fem.variables.d_u;

%% P1 to P0
aa.connec = s.mesh.connec;
aa.type   = s.mesh.type;
aa.fNodes = u;
fefunDisp = FeFunction(aa);
p1displac = fefunDisp.computeValueInCenterElement();

%% P0 to P1
bb.mesh   = s.mesh;
bb.connec = s.mesh.connec;
bb.nelem  = size(s.mesh.connec,1);
bb.nnode  = size(s.mesh.connec,2);
bb.npnod  = size(s.mesh.coord,1);
projector = Projector_P0toP1(bb);

strainP1 = projector.project(strain);