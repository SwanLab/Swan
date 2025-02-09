% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique
% Algiane Froehly, CARDAMOME, INRIA-SOFT
% LGPL Lesser General Public License v3.0.
% Remeshing using Mmg tools : https://www.mmgtools.org

clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N        = 1e2;
filename = 'testIO.mesh';
Ireq     = [2 47]';
h        = 0.1;


%%% TRIANGLES
% Generate surfacic mesh
tic
mesh1 = mshDisk(N,1);
toc

% Add colours
ctr                    = mesh1.ctr;
mesh1.col(ctr(:,1)>=0) = 1;
mesh1.col(ctr(:,1)<0)  = 2;

% Graphical representation
figure
plot(mesh1)
hold on
plot(mesh1.sub(Ireq),'r')
axis equal
colorbar

% Write mesh file
tic
mmgMeshWrite(filename,mesh1,Ireq)
toc

% Read mesh file
tic
mesh2 = mmgMeshRead(filename);
toc

% Validation 
norm(mesh1.vtx-mesh2.vtx,'inf')
norm(mesh1.elt-mesh2.elt,'inf')
norm(mesh1.col-mesh2.col,'inf')

% Mmg refinment
refine = mmg(mesh1);
req(refine,Ireq);
map(refine,h*ones(size(mesh1.vtx,1),1))
verbose(refine,0)
tic
mesh3 = run(refine);
toc

% Intersection 
inter = intersect(mesh3,mesh1)

% Graphical representation
figure
plot(mesh3)
hold on
plot(inter,'r')
axis equal
colorbar


%%% TETRA
% Generate surfacic mesh
tic
mesh1 = mshCube(N,[1 1 1]);
toc

% Add colours
ctr                    = mesh1.ctr;
mesh1.col(ctr(:,1)>=0) = 1;
mesh1.col(ctr(:,1)<0)  = 2;

% Graphical representation
figure
plot(mesh1)
hold on
plot(mesh1.sub(Ireq),'r')
alpha(0.1)
axis equal
colorbar

% Write mesh file
tic
mmgMeshWrite(filename,mesh1,Ireq)
toc

% Read mesh file
tic
mesh2 = mmgMeshRead(filename);
toc

% Validation 
norm(mesh1.vtx-mesh2.vtx,'inf')
norm(mesh1.elt-mesh2.elt,'inf')
norm(mesh1.col-mesh2.col,'inf')

% Mmg refinment
refine = mmg(mesh1);
req(refine,Ireq);
map(refine,h*ones(size(mesh1.vtx,1),1))
verbose(refine,0)
tic
mesh3 = run(refine);
toc

% Intersection
inter = intersect(mesh3,mesh1)

% Graphical representation
figure
plot(mesh3)
hold on
plot(inter,'r')
alpha(0.1)
axis equal
colorbar






disp('~~> Michto gypsilab !')



