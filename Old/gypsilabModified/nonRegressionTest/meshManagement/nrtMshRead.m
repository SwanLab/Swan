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
%|    #    |   FILE       : nrtMshRead.m                                  |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Mesh generation and readers                   |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Segment
Nvtx = 1e1;
L    = 1;
mesh = mshSegment(Nvtx,L);
figure
plot(mesh)
axis equal
view(30,30)

% Square
Nvtx = 1e2;
L    = [1 1];
mesh = mshSquare(Nvtx,L);
figure
plot(mesh)
axis equal
view(30,30)

% Disk
Nvtx = 1e2;
rho  = 2;
mesh = mshDisk(Nvtx,rho);
figure
plot(mesh)
axis equal
view(30,30)

% Cube
Nvtx = 1e3;
L    = [1 1 1];
mesh = mshCube(Nvtx,L);
figure
plot(mesh)
axis equal
view(30,30)
alpha(0.3)

% Sphere
Nvtx = 1e3;
rho  = 2;
mesh = mshSphere(Nvtx,rho);
figure
plot(mesh)
axis equal
view(30,30)
alpha(0.3)

% .msh triangle
mesh = msh('sphere.msh');
figure
plot(mesh)
axis equal
view(30,30)
alpha(0.3)

% .msh tetra
mesh = msh('sphereTet.msh');
figure
plot(mesh)
axis equal
view(30,30)
alpha(0.3)

% .ply
mesh = msh('cube.ply');
figure
plot(mesh)
axis equal
view(30,30)
alpha(0.3)

% .stl
mesh = msh('cube.stl');
figure
plot(mesh)
axis equal
view(30,30)
alpha(0.3)

% .mesh
mesh = msh('cube.mesh');
figure
plot(mesh)
axis equal
view(30,30)
colorbar
alpha(0.3)

% .vtk
mesh = msh('plan.vtk');
figure
plot(mesh)
view(90,0)



disp('~~> Michto gypsilab !')




