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
%|    #    |   FILE       : nrtMshClean.m                                 |
%|    #    |   VERSION    : 0.53                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2019                                    |
%| ( === ) |   SYNOPSIS   : Clean tetrahedral mesh                        |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Load mesh
load tetmesh

% Build mesh
mesh1 = msh(X,tet);

% Graphical representation
figure
plot(mesh1)
axis equal
view(45,45)

% Add vertices to tetrahedral mesh
tmp = [zeros(10,3) ; 10*rand(10,3)];
tet = tet + size(tmp,1);
X   = [tmp ; X ; 10*rand(20,3)];

% Build and clean mesh
mesh2 = msh(X,tet);

% Test mesh egality
~isequal(mesh1,mesh2)

% Graphical representation
figure
plot(mesh2)
axis equal
view(45,45)

% Colours
mesh = mshSquare(20,[1 1]);
ctr  = mesh.ctr;
mesh.col(ctr(:,1)<0)  = 1;    
mesh.col(ctr(:,1)>=0) = 2;

% Graphical representation
figure
plot(mesh)
axis equal
view(45,45)

% Sub-meshing with small translation 
delta     = 1e-2;
mesh1     = mesh.sub(mesh.col==1);
mesh2     = mesh.sub(mesh.col==2);
mesh2.vtx = mesh2.vtx + delta;
meshT     = union(mesh1,mesh2);

% Graphical representation
figure
plot(meshT)
axis equal
view(45,45)

% Clean with specified range
stp   = mesh.stp;
meshT = clean(meshT,0.5*stp(1));




disp('~~> Michto gypsilab !')


