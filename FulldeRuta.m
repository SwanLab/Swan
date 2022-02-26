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
x = 0:0.5:1;
y = 0:0.5:1;
z = 0:0.5:1;
[X,Y,Z] = meshgrid(x,y,z);
% X = reshape(X, 27,1);
% Y = reshape(Y, 27,1);
% Z = reshape(Z, 27,1);

% plot3(X(:),Y(:),Z(:),'ok')

fvc = surf2patch(reshape(X,9,3),reshape(Y,9,3),reshape(Z,9,3));
% fvc.vertices(:,3) = []; % 2D

coords = fvc.vertices;
connec = fvc.faces;
%% 3D Example
% Cantilever beam
clc; clear; close all;

% Coordinates
x = 0:0.1:1;
y = 0:0.05:0.25;
z = 0:0.05:0.25;

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

%% Todo
% {{done}} Create example 2D not using load
% Create example 3D not using load
% eliminate istre,jstre loop 
% investigate how to efficiently multiply B,C,B
% {{done}} Use BmatrixComputer in LHSintegrator_StifnessElastic
% With large example compare Sparse vs Accumarray
% eliminate computeLHS from Integrator_Simple
% Element_DiffReact K, M, Mr with LHSintegrator
% Eliminate computeLHS from integratorComposite

%% Endgame
% 7. Force use integrator for assembly
% 8. LHS integrator must be composed by integrator_simple for assemnbly