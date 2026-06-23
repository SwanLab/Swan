clear;
close all;

%% COMPOSITION OF LEVELSETS 2D
x1       = linspace(0,6,40);
x2       = linspace(0,6,40);
[xv,yv]  = meshgrid(x1,x2);
[F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
m.coord  = V(:,1:2);
m.connec = F;
mesh     = Mesh.create(m);
bmesh= mesh.createBoundaryMesh();

r=[  0.3,0.2,0.5
     0.4,0.7,0.2
     0.6,0.9,0.8];

length=[2,1,0.5
        1.5,0.9,0.2
        2,1.3,0.5];

[x0,y0]=ComputeSubdomainCentroids(r,mesh);

[Nx, Ny] = size(r);

% type =repmat("Circle",Nx,Ny);
% 
% type=["Circle","Square","Circle"
%       "CrossedSquare","Circle","Circle"
%       "Circle","Square","Square"];

type=["CrossedSquare","Square","Circle"
      "Square","Square","Circle"
      "CrossedSquare","Square","Square"];
paramMatrix(Nx,Ny) = struct('type',[],'radius',[],'xCoorCenter',[],'yCoorCenter',[],'length',[],'tFrame',[],...
                            'tCross',[]);

for i = 1:Nx
    for j = 1:Ny
        paramMatrix(i,j).type        = type(i,j);
        paramMatrix(i,j).type        = "CrossedSquare";
        paramMatrix(i,j).radius      = r(i,j);
        paramMatrix(i,j).length      = 2;
        paramMatrix(i,j).tCross      = 0.2;
        paramMatrix(i,j).tFrame      = 0.2;
        paramMatrix(i,j).xCoorCenter = x0(i,j);
        paramMatrix(i,j).yCoorCenter = y0(i,j);
    end
end

% paramMatrix(1,1).tCross=0.2;
% paramMatrix(1,1).tFrame=0.2;
% paramMatrix(3,1).tCross=0.3;
% paramMatrix(3,1).tFrame=0.3;

% params.r  = @(x) paramFromGrid(x,r,mesh);
% params.x0 = @(x) paramFromGrid(x,x0,mesh);
% params.y0 = @(x) paramFromGrid(x,y0,mesh);

s.type       = 'GivenPattern';
s.paramsList = paramMatrix;
g        = GeometricalFunction(s);
phiFun   = g.computeLevelSetFunction(mesh);
obj.levelSet = phiFun;
ls       = phiFun.fValues;


%% COMPOSITION OF LEVELSETS 3D

mesh= TetraMesh(6,6,6,20,20,20);
r=zeros(3,3,3);
r(:,:,1)=[  0.2,0.5,0.8
            0.2,0.5,0.8
            0.2,0.5,0.8];
r(:,:,2)=r(:,:,1);
r(:,:,3)=r(:,:,1);

[Nx, Ny, Nz] = size(r);

[x0,y0,z0]=ComputeSubdomainCentroids3D(r,mesh);

paramMatrix(Nx,Ny,Nz) = struct('type',[],'radius',[],'xCoorCenter',[],'yCoorCenter',[],'zCoorCenter',[]);

for i = 1:Nx
    for j = 1:Ny
        paramMatrix(i,j).type        = type(i,j);
        paramMatrix(i,j).type        = "CrossedSquare";
        paramMatrix(i,j).radius      = r(i,j);
        paramMatrix(i,j).length      = 2;
        paramMatrix(i,j).tCross      = 0.2;
        paramMatrix(i,j).tFrame      = 0.2;
        paramMatrix(i,j).xCoorCenter = x0(i,j);
        paramMatrix(i,j).yCoorCenter = y0(i,j);
    end
end

%% REFERENCE LEVEL SET
x1      = linspace(-1,1,20);
x2      = linspace(-1,1,20);
[xv,yv] = meshgrid(x1,x2);
[F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
m.coord  = V(:,1:2);
m.connec = F;
mesh     = Mesh.create(m);

% gPar.type         = 'LatticeCircle';
% gPar.length       = 2;
% gPar.xCoorCenter  = 0;
% gPar.yCoorCenter  = 0;
% gPar.tFrame       = 0.15;
% gPar.tCross       = 0.15;
% gPar.meanRadius   = 0.5;
% gPar.tRadius      = 0.15;
% 
% 
% gPar.type         = 'LatticeCircleV2';
% gPar.length       = 2;
% gPar.xCoorCenter  = 0;
% gPar.yCoorCenter  = 0;
% gPar.tFrame       = 0.15;
% gPar.tCross       = 0.15;
% gPar.radius       = 2;
% 
% gPar.type       = 'Square';
% gPar.length     =   0.5 ;
% gPar.xCoorCenter  =  0.5;
% gPar.yCoorCenter  =  0.5;
% 
% gPar.type         = 'CircleInclusion';
% gPar.radius       = 0.5;
% gPar.xCoorCenter  = 0;
% gPar.yCoorCenter  = 0;


g        = GeometricalFunction(gPar);
phiFun   = g.computeLevelSetFunction(mesh);
obj.levelSet = phiFun;
ls       = phiFun.fValues;

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(ls);
uMesh.plot();



%% AUXETIC CELL

x1      = linspace(-1.5,1.5,50);
x2      = linspace(-1,1,50);
[xv,yv] = meshgrid(x1,x2);
[F,V]   = mesh2tri(xv,yv,zeros(size(xv)),'x');
m.coord  = V(:,1:2);
m.connec = F;
mesh     = Mesh.create(m);

gPar.type           = 'Auxetic';
gPar.length         = 2;
gPar.height         = 2;
gPar.xCoorCenter    = 0;
gPar.yCoorCenter    = 0;
gPar.theta          = 60;
gPar.beta           = 64;
gPar.thickness      = 0.2;


g        = GeometricalFunction(gPar);
phiFun   = g.computeLevelSetFunction(mesh);
obj.levelSet = phiFun;
ls       = phiFun.fValues;

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(ls);
uMesh.plot();

%% lattice 3D
mesh=TetraMesh(1,1,1,40,40,40);

gPar.type         = 'CrossedSquare3D';
gPar.length       = 2;
gPar.xCoorCenter  = 0;
gPar.yCoorCenter  = 0;
gPar.zCoorCenter  = 0;
gPar.yCoorCenter  = 0;
gPar.tFrame       = 0.2;
gPar.tCross       = 0.2;

g        = GeometricalFunction(gPar);
phiFun   = g.computeLevelSetFunction(mesh);
obj.levelSet = phiFun;
ls       = phiFun.fValues;

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(-ls);
uMesh.plot();

mS=uMesh.createInnerMesh();
mS.print("LatticeMesh","paraview");


%% Cube 3D

mesh=TetraMesh(1,1,1,10,10,10);

gPar.type         = 'Cube';
gPar.length       = 2;
gPar.xCoorCenter  = 0;
gPar.yCoorCenter  = 0;
gPar.height       = 1;


g        = GeometricalFunction(gPar);
phiFun   = g.computeLevelSetFunction(mesh);
obj.levelSet = phiFun;
ls       = phiFun.fValues;

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(ls);
uMesh.plot();

mS=uMesh.createInnerMesh();
mS.print("LatticeMesh","paraview");

%% NACA 3D

lOld= 8;
hOld=4;
nxOld=420;
nyOld=round(nxOld / lOld * hOld / 0.8);
hy = hOld/nyOld;
aspect=lOld/hOld;

Mesh.length   = 1.05;
Mesh.height   = 0.2;
scale         = lOld/Mesh.length;
Mesh.nx       = round(nxOld/scale);
Mesh.ny       = nyOld * (Mesh.height / hOld);
% refMesh = TriangleMesh(Mesh.length,Mesh.height,Mesh.nx,Mesh.ny);
% refMesh = TriangleMesh(1.05,Mesh.height,60,11);


refMesh= TetraMesh(Mesh.length,Mesh.height,0.5,60,20,10);
gPar.type     = 'Naca';
gPar.m        = 0.022;
gPar.p        = 0.4;
gPar.t        = 0.15;
gPar.chord    = 1;
gPar.AoA      = 5;
gPar.xLE      = (Mesh.length - gPar.chord) / 2;
gPar.yLE      = Mesh.height/2;
gPar.wall     = 0.02;
g        = GeometricalFunction(gPar);
phiFun   = g.computeLevelSetFunction(refMesh);
obj.levelSet = phiFun;
ls       = phiFun.fValues;

sUm.backgroundMesh = refMesh;
sUm.boundaryMesh   = refMesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(ls);
uMesh.plot();
% 
mS=uMesh.createInnerMesh();
mS.print("AirfoilMesh","paraview");
%% Functions
function [x0,y0] = ComputeSubdomainCentroids(param,mesh)
    Nx = size(param,2);
    Ny = size(param,1);
    xmin=min(mesh.coord(:,1));
    xmax=max(mesh.coord(:,1));
    ymin=min(mesh.coord(:,2));
    ymax=max(mesh.coord(:,2));
    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;

    x_center = xmin + dx/2 : dx : xmax - dx/2;
    y_center = ymin + dy/2 : dy : ymax - dy/2;

    [x0, y0] = meshgrid(x_center, y_center);
end

function [x0,y0,z0] = ComputeSubdomainCentroids3D(param,mesh)
    Nx = size(param,2);
    Ny = size(param,1);
    Nz = size(param,3);
    xmin = min(mesh.coord(:,1));
    xmax = max(mesh.coord(:,1));
    ymin = min(mesh.coord(:,2));
    ymax = max(mesh.coord(:,2));
    zmin = min(mesh.coord(:,3));
    zmax = max(mesh.coord(:,3));
    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;
    dz = (zmax - zmin)/Nz;
    x_center = xmin + dx/2 : dx : xmax - dx/2;
    y_center = ymin + dy/2 : dy : ymax - dy/2;
    z_center = zmin + dz/2 : dz : zmax - dz/2;
    [x0, y0, z0] = ndgrid(x_center, y_center, z_center);

end