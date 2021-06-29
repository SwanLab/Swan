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
%|    #    |   FILE       : nrtMshBitree.m                                |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Binary tree subdivision for tetra, triangles, |
%|  `---'  |                edges and particles                           |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Initialize 3D mesh with tetra
Nvtx  = 1e3;
L     = [4 3 2];
tet   = mshCube(Nvtx,L)
Nleaf = 1;

% Tetra tree
tree(tet,'binary',Nleaf,2);
tic
tree(tet,'binary');
toc

% Triangle tree
tri = tet.bnd
tree(tri,'binary',Nleaf,3);
tic
tree(tri,'binary');
toc

% Edge tree
edg = tri.edg
tree(edg,'binary',Nleaf,4);
tic
tree(edg,'binary');
toc

% Particles tree
prt = tri.prt
tmp = tree(prt,'binary',Nleaf,5);
tic
tree(prt,'binary');
toc

% Error
ind      = cell2mat(tmp{end}.ind);
X        = tmp{end}.ctr;
X(ind,:) = X;
norm(prt.vtx-X)

% Non uniform mesh
mesh = mshSphere(1e3,1);
fct  = @(X) floor(3*(1+X(:,1))/2);
refi = mesh.refine(fct)
tree(refi,'binary',Nleaf,6);
tic
tree(refi,'binary');
toc


disp('~~> Michto gypsilab !')


