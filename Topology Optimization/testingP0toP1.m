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
z.fValues = strain(:,1:2);
z.connec  = s.mesh.connec;
z.type    = s.mesh.type;
strainFun = P0Function(z);

strainP1Disc = strainFun.computeP1DiscontinuousFunction();
% strainP1Disc.plot(s.mesh);

%% Create projector
bb.mesh   = s.mesh;
bb.connec = s.mesh.connec;
projector = Projector_P0toP1(bb);
p1strain = projector.project(strainFun);

%% toP1disc
projectorDisc = ProjectorToP1discont(bb);
ss.origin = 'P0';
ss.x      = strainFun;
strainFundisc = projectorDisc.project(ss);

%% Lets create a FGaussFun
x.fValues = permute(fem.variables.strain,[2 1 3]);
x.connec  = s.mesh.connec;
x.type    = s.mesh.type;
strainGaussDiscFun = FGaussDiscontinuousFunction(x);

%% Use this to create results for testing:
nrows = numel(p1strain.fValues);
xP = reshape(p1strain.fValues,[nrows,1]);
% save('name.mat','xP');