%% Filter

clear;
clc;
close all;

% filename   = 'anisoCantilever';
% a.fileName = filename;
% gid        = FemDataContainer(a);
% mesh       = gid.mesh;
% h          = mesh.computeMinCellSize();
% 
% s.fHandle = @(x) 1-heaviside((x(1,:,:)-1).^2+(x(2,:,:)-0.5).^2-0.3.^2);
% s.ndimf   = 1;
% s.mesh    = mesh;
% fun{1}    = AnalyticalFunction(s);

mesh = TriangleMesh(1,1,5,5);

s.fHandle = @(x) 1-heaviside((x(1,:,:)-0.5).^2+(x(2,:,:)-0.5).^2-0.3.^2);
s.ndimf   = 1;
s.mesh    = mesh;
fun{1}    = AnalyticalFunction(s);

s.mesh  = mesh;
s.alpha = 4;
s.beta  = 0;
s.theta = 90;
filter  = NonLinearFilterSegment(s);
filter.updateEpsilon(1);

[fun{end+1},err] = filter.compute(fun{1},2);

figure
plot(err)
fun{end}.plot();


%% Cantilever

clear;
clc;
close all;

load('Experiments/Density/CantileverDensityOriginalDesVar.mat');
mesh = TriangleMesh(2,1,100,50);
h    = mesh.computeMinCellSize();


s.mesh  = mesh;
s.alpha = 4;
s.beta  = 0;
s.theta = 90;
s.tol0  = 1e-6;
filter  = NonLinearFilterSegment(s);

levelSet         = -ones(mesh.nnodes,1);
sG.backgroundMesh = mesh;
sG.boundaryMesh   = mesh.createBoundaryMesh();
uMesh = UnfittedMesh(sG);
uMesh.compute(levelSet);

ss.mesh    = mesh;
ss.uMesh   = uMesh;
ss.filter  = filter;
ss.epsilon = h;
ss.value0  = 6;
pF         = PerimeterFunctional(ss);

sD.plotting = false;
sD.type = 'Density';
sD.fun = f;
x = Density(sD);
p          = pF.computeFunctionAndGradient(x); % 4.1