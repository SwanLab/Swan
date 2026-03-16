clear;
close all;

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

% gPar.type         = 'CrossedSquare';
% gPar.length       = 1;
% gPar.xCoorCenter  = 0.5;
% gPar.yCoorCenter  = 0.5;
% gPar.tFrame       = 0.03;
% gPar.tCross       = 0.7;
% 
% gPar.type       = 'Square';
% gPar.length     =   0.5 ;
% gPar.xCoorCenter  =  0.5;
% gPar.yCoorCenter  =  0.5;
% 
% gPar.type         = 'Circle';
% gPar.radius       = params.r;
% gPar.xCoorCenter  = params.x0;
% gPar.yCoorCenter  = params.y0;

% g        = GeometricalFunction(gPar);

% s.fBase  = g.getHandle;
s.type       = 'GivenPattern';
s.paramsList = paramMatrix;
% s.r     = r;
% s.x0    =x0;
% s.y0    =y0;
g        = GeometricalFunction(s);
phiFun   = g.computeLevelSetFunction(mesh);
obj.levelSet = phiFun;
ls       = phiFun.fValues;



sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(-ls);
uMesh.plot();


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