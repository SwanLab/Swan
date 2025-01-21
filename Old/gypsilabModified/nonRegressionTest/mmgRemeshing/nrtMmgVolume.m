% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique  
% Algiane Froehly, CARDAMOME, INRIA-SOFT 
% LGPL Lesser General Public License v3.0. 
% Remeshing using Mmg tools : https://www.mmgtools.org         

% Cleaning
clear all
close all
clc

% Library path
run('../../addpathGypsilab')

% Mesh unit square
mesh = msh('torusholes.msh');

% Colours
ctr      = mesh.ctr;     
mesh.col = ceil(3+1/max(ctr(:,1)) .* ctr(:,1)); 

% Distance map
stp = mesh.stp;
I   = abs(mesh.vtx(:,1));
I   = 10*stp(3)/max(abs(I)) .* I;

% Anisotropic
Z = zeros(size(I));
U = ones(size(I));
I = 1./(I.^2);
I = [I,Z,I,Z,Z,I];

% Graphical representation
figure
plot(mesh)
axis equal
colorbar
xlabel('X')

% Refinment with Haussdorf
refine  = mmg(mesh,1e-2);
[mesh1] = run(refine);

% Graphical representation
figure
plot(mesh1)
hold on
axis equal
colorbar

% Options
% aniso(refine)
% refine.aniso()
% hgrad(refine,-1)
% hmin(refine,0.1)
% hmax(refine,0.2)
% hsiz(refine,0.05)
% hausd(refine,1e-2) 
% noinsert(refine)
% nomove(refine)
% noswap(refine)
% angle(refine,90)
% memory(refine,1e3)
verbose(refine,0)
% options(refine)
% refine.mesh(mesh)
map(refine,I)

% options(refine)
mesh2 = run(refine);

% Graphical representation
figure
plot(mesh2)
hold on
axis equal
colorbar


    

disp('~~> Michto gypsilab !')




