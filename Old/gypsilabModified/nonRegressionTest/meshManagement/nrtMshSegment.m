%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtMshSegment.m                               |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Edge mesh of a segment                        |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Create mesh
Nvtx = 21;
L    = 4;
mesh = mshSegment(Nvtx,L);

% Colours
ctr = 1/2 .* (...
        mesh.vtx(mesh.elt(:,1),:) + ...
        mesh.vtx(mesh.elt(:,2),:) ) ;
mesh.col(ctr(:,1)<0)               = 1;    
mesh.col(ctr(:,1)>=0 & ctr(:,1)<1) = 2;
mesh.col(ctr(:,1)>=1)              = 3;

% Graphical representation
figure
plot(mesh)
colorbar

% Sub-meshing
mesh1 = mesh.sub(mesh.col==1);
mesh2 = mesh.sub(mesh.col==2);
mesh3 = mesh.sub(mesh.col==3);
norm(mesh1.col-1,'inf')
norm(mesh2.col-2,'inf')
norm(mesh3.col-3,'inf')

% Graphical representation
figure
plot(mesh1)
hold on
plot(mesh2)
plot(mesh3)
colorbar

% Center
Xctr = mesh.ctr;
hold on
plot3(Xctr(:,1),Xctr(:,2),Xctr(:,3),'*y')

% Normals
plotNrm(mesh,'r')

% Volume
V = mesh.ndv;
norm(sum(V)-L)/L

% Length
l = mesh.stp;

% Edges
figure
plot(mesh.edg)
colorbar

% Particles
figure
plot(mesh.prt)
colorbar

% Unicity
tmp     = mesh;
tmp.vtx = [tmp.vtx ; tmp.vtx];
tmp.elt = [tmp.elt ; tmp.elt + size(mesh.vtx,1)];
tmp.col = [tmp.col ; tmp.col];
[~,I]   = unique(tmp);
norm(I-(1:size(mesh.elt,1))','inf')

% Intersection
[~,I] = intersect(mesh1,mesh);
norm(I-(1:size(mesh1.elt,1))','inf')

% Union
[tmp,I] = union(mesh,mesh1);
norm(I-(1:size(mesh.elt,1))','inf')

% Difference
[mesh4,~] = setdiff(mesh,mesh1);
figure
plot(mesh4)
colorbar

% Reconstruction
mesh5 = union(mesh1,mesh2);
mesh5 = union(mesh5,mesh3);
setdiff(mesh,mesh5)

% Boundary
bound = mesh5.bnd;
figure
plot(bound)
colorbar

% Clean degenerated mesh
mesh6 = mesh;
ind   = (mesh6.vtx(:,1)>0);
mesh6.vtx(ind,1) = 0;
mesh6 = clean(mesh6,1e-6);
setdiff(mesh6,mesh1)

% Clean degenerated particle mesh
vtx = mesh.vtx;
vtx(vtx(:,1)>0,1) = 0;
mesh7 = msh(vtx);
mesh7 = clean(mesh7,1e-6);
setdiff(mesh7,msh(mesh1.vtx))



disp('~~> Michto gypsilab !')
