%% 2D Example
% Cantilever beam
clc; clear; close all;

% Coordinates
x = 0:0.1:1;
y = 0:0.05:0.25;

%Mesh
[X,Y,Z] = meshgrid(x,y,0);
fvc = surf2patch(X,Y,Z,'triangles');
fvc.vertices(:,3) = []; % 2D
coords = fvc.vertices;
connec = fvc.faces;
m.coord = coords;
m.connec = connec;
p.mesh = Mesh(m);
npnod = size(coords,1);

% Boundary conditions
dirichNodes = find(coords(:,1) == 0 );
ndirich = size(dirichNodes,1);
dirichlet = [dirichNodes,   ones(ndirich,1), zeros(ndirich,1);
             dirichNodes, 2*ones(ndirich,1), zeros(ndirich,1)];
neumann   = [npnod, 2, 0.0001];

% FEM parameters
p.dim = '2D';
p.type = 'ELASTIC';
p.scale = 'MACRO';
p.dirichlet = dirichlet;
p.pointload = neumann;

% Solution
fem = NewFEM.create(p);
fem.solve();
fem.plot();

%% 3D reference
clc; clear; close all;
load("newFemHexahedra.mat")

%% Small 3D attempt
clc; clear;

cm.coord = [0, 0, 0;
            1, 0, 0;
            0, 1, 0;
            1, 1, 0;
            0, 0, 1;
            1, 0, 1;
            0, 1, 1;
            1, 1, 1;
           ];
cm.connec = [1 2 3 4 5 6 7 8];
cubeMesh = Mesh(cm);
dirichlet = [1 1 0;
             1 2 0;
             1 3 0;
             2 1 0;
             2 2 0;
             2 3 0;
             5 1 0;
             5 2 0;
             5 3 0;
             6 1 0;
             6 2 0;
             6 3 0;];
neumann = [7, 3, 0.5;
           8, 3, 0.5;
           ];
p.dim = '3D';
p.type = 'ELASTIC';
p.scale = 'MACRO';
p.mesh  = cubeMesh;
p.dirichlet = dirichlet;
p.pointload = neumann;

fem = NewFEM.create(p);
fem.solve();
fem.plot();

%% 3D Example
% Cantilever beam
clc; clear; close all;

% Coordinates
x = 0:0.1:1;
y = 0:0.05:0.25;
z = 0:0.05:0.25;

% x = linspace(0,1,10);
% y = linspace(0,1,10);
% z = linspace(0,2,20);
% [X,Y,Z] = meshgrid(x,y,z);
% coord  = [X(:) Y(:) Z(:)];
% d = delaunayTriangulation(coord);
% s.connec = d.ConnectivityList;
% s.coord  = coord;
% obj.backgroundMesh = Mesh(s);

%Mesh
[X,Y,Z] = meshgrid(x,y,z);
% fvc = surf2patch(X,Y,Z,'f');
sizeX = size(X,1);
sizeY = size(X,2);
sizeZ = size(X,3);
fvc = surf2patch(reshape(X,sizeX*sizeY,sizeZ),reshape(Y,sizeX*sizeY,sizeZ),reshape(Z,sizeX*sizeY,sizeZ));
% fvc.vertices(:,3) = []; % 2D

coords = fvc.vertices;
connec = fvc.faces;
m.coord = coords;
m.connec = connec;
p.mesh = Mesh(m);
npnod = size(coords,1);

% Boundary conditions
dirichNodes = find(coords(:,1) == 0 );
ndirich = size(dirichNodes,1);
dirichlet = [dirichNodes,   ones(ndirich,1), zeros(ndirich,1);
             dirichNodes, 2*ones(ndirich,1), zeros(ndirich,1);
             dirichNodes, 3*ones(ndirich,1), zeros(ndirich,1)];
neumann   = [npnod, 2, 1];

% FEM parameters
p.dim = '3D';
p.type = 'ELASTIC';
p.scale = 'MACRO';
p.dirichlet = dirichlet;
p.pointload = neumann;

% Solution
fem = NewFEM.create(p);
fem.solve();
fem.plot();

%% Visualize

%             Tn = connec;
%             x  = coords(:,1);
%             y  = coords(:,2);
%             figure()
%             hold on
%             colormap jet;
%             plot(x(Tn)',y(Tn)','--','linewidth',0.5);
