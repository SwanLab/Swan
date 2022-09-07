%% Create sample FEM results
clear; close all;

file = 'test2d_triangle';
a.fileName = file;
s = FemDataContainer(a);
fem = FEM.create(s);
fem.solve();

strain = fem.variables.strain;
uCol   = fem.variables.d_u;
u      = reshape(uCol,[s.mesh.ndim,s.mesh.nnodes])';

%% P1 to P0
aa.connec = s.mesh.connec;
aa.type   = s.mesh.type;
aa.fNodes = u;
fefunDisp = FeFunction(aa);
p0displac = fefunDisp.computeValueInCenterElement()';

%% P0 to P1
bb.mesh   = s.mesh;
bb.connec = s.mesh.connec;
bb.nelem  = size(s.mesh.connec,1);
bb.nnode  = size(s.mesh.connec,2);
bb.npnod  = size(s.mesh.coord,1);
projector = Projector_P0toP1(bb);

% strainXX = squeeze(strain(:,1,:));
% strainP1 = projector.project(strainXX);

% strainOnes = ones(16,1);
% strainOnesP1 = projector.project(strainOnes);

strainCol = reshape(strain, [bb.nelem*3, 1]);
strainP1 = projector.project(strainCol);

%% P1 to P0 v2
cc.mesh   = s.mesh;
cc.connec = s.mesh.connec;
cc.nelem  = size(s.mesh.connec,1);
cc.nnode  = size(s.mesh.connec,2);
cc.npnod  = size(s.mesh.coord,1);
projector2 = Projector_P1toP0(cc);
u_P0 = projector2.project(u);
error = max(max(abs((p0displac-u_P0)./p0displac)));