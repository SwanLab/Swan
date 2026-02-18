clear;
close all;

x1       = linspace(0,1,40);
x2       = linspace(0,1,40);
[xv,yv]  = meshgrid(x1,x2);
[F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
m.coord  = V(:,1:2);
m.connec = F;
mesh     = Mesh.create(m);

gPar.type         = 'CrossedSquare';
gPar.length       = 1;
gPar.xCoorCenter  = 0.5;
gPar.yCoorCenter  = 0.5;
gPar.tFrame       = 0.1;
gPar.tCross       = 0.1;
% 
% gPar.type       = 'Square';
% gPar.length     =   0.5 ;
% gPar.xCoorCenter  =  0.5;
% gPar.yCoorCenter  =  0.5;

% gPar.type         = 'Circle';
% gPar.radius       = 0.25;
% gPar.xCoorCenter  = 0.5;
% gPar.yCoorCenter  = 0.5;

g                 = GeometricalFunction(gPar);
phiFun            = g.computeLevelSetFunction(mesh);
lsCircle          = phiFun.fValues;
lsCircleInclusion = -lsCircle;

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(lsCircleInclusion);

uMesh.plot;