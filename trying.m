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
ss.epsilon = 2*h;
ss.value0  = 1;
pF         = PerimeterFunctional(ss);

sD.plotting = false;
sD.type = 'Density';
sD.fun = f;
x = Density(sD);
p          = pF.computeFunctionAndGradient(x); % 17.2675


%% Iso and Ani

clear;
clc;
close all;

mesh = TriangleMesh(1,1,200,200);
h    = mesh.computeMinCellSize();

s.type = 'CircleInclusion';
s.length = 0.5;
s.xCoorCenter = 0.5;
s.yCoorCenter = 0.5;
s.radius = 0.25;
g = GeometricalFunction(s);
ls = g.computeLevelSetFunction(mesh);

sU.backgroundMesh = mesh;
sU.boundaryMesh = mesh.createBoundaryMesh();
uMesh = UnfittedMesh(sU);
uMesh.compute(ls.fValues);

chi = CharacteristicFunction.create(uMesh);

u = 80;
alpha = 90;
CAnisotropic = [tand(u),0;0,1/tand(u)];
R = [cosd(alpha),-sind(alpha)
    sind(alpha), cosd(alpha)];
CGlobal = R*CAnisotropic*R';

sF.A     = ConstantFunction.create(CGlobal,mesh);
sF.trial = LagrangianFunction.create(mesh,1,'P1');
sF.mesh = mesh;
sF.filterType = 'PDE';
sF.boundaryType = 'Neumann';
sF.metric = 'Isotropy';
filter = Filter.create(sF);
filter.updateEpsilon(10*h);

rhoEps = filter.compute(chi,3);
rhoEps.print('RhoEps')

fun = chi.*(1-rhoEps);

sFF.trial = LagrangianFunction.create(mesh,1,'P1');
sFF.mesh = mesh;
filter = FilterLump(sFF);
int = filter.compute(fun,3);
int.print('Integrand');


%% Segment

clear;
clc;
close all;

mesh = TriangleMesh(1,1,200,200);
h    = mesh.computeMinCellSize();

s.type = 'SquareInclusion';
s.length = 0.5;
s.xCoorCenter = 0.5;
s.yCoorCenter = 0.5;
s.radius = 0.25;
g = GeometricalFunction(s);
ls = g.computeLevelSetFunction(mesh);

sU.backgroundMesh = mesh;
sU.boundaryMesh = mesh.createBoundaryMesh();
uMesh = UnfittedMesh(sU);
uMesh.compute(ls.fValues);

chi = CharacteristicFunction.create(uMesh);

sF.mesh  = mesh;
sF.alpha = 4;
sF.beta  = 0;
sF.theta = 90;
sF.tol0  = 1e-6;
filter  = NonLinearFilterSegment(sF);
filter.updateEpsilon(10*h);

rhoEps = filter.compute(chi,3);
rhoEps.print('RhoEps')

fun = chi.*(1-rhoEps);

sFF.trial = LagrangianFunction.create(mesh,1,'P1');
sFF.mesh = mesh;
filter = FilterLump(sFF);
int = filter.compute(fun,3);
int.print('Integrand');