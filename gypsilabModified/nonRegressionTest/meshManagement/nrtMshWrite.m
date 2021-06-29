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
%|    #    |   FILE       : nrtMshWrite.m                                 |
%|    #    |   VERSION    : 0.53                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : Mesh readers and writer                       |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')


%%% Particles mesh
mesh = msh(rand(100,3));

% Graphical representation
figure
plot(mesh)
axis equal
view(10,10)
grid on

% .msh format
mshWriteMsh('testPoint.msh',mesh)
meshMsh = msh('testPoint.msh');
norm(meshMsh.vtx - mesh.vtx)
norm(meshMsh.elt - mesh.elt)


%%% Edge mesh
mesh = mshSegment(100,1);

% Graphical representation
figure
plot(mesh)
axis equal
view(10,10)
grid on

% .msh format
mshWriteMsh('testEdge.msh',mesh)
meshMsh = msh('testEdge.msh');
norm(meshMsh.vtx - mesh.vtx)
norm(meshMsh.elt - mesh.elt)


%%% Surfacic mesh
mesh = mshSphere(100,1);

% Graphical representation
figure
plot(mesh)
axis equal
view(30,30)

% .ply format
mshWritePly('testTriangle.ply',mesh)
meshPly = msh('testTriangle.ply');
norm(meshPly.vtx - mesh.vtx)
norm(meshPly.elt - mesh.elt)

% .vtk format
mshWriteVtk('testTriangle.vtk',mesh)
meshVtk = msh('testTriangle.vtk');
norm(meshVtk.vtx - mesh.vtx)
norm(meshVtk.elt - mesh.elt)

% .msh format
mshWriteMsh('testTriangle.msh',mesh)
meshMsh = msh('testTriangle.msh');
norm(meshMsh.vtx - mesh.vtx)
norm(meshMsh.elt - mesh.elt)


%%% Volumic mesh
mesh = mshCube(100,[1 1 1]);

% Graphical representation
figure
plot(mesh)
axis equal
view(10,10)
alpha(0.3)

% .msh format
mshWriteMsh('testTetra.msh',mesh)
meshMsh = msh('testTetra.msh');
norm(meshMsh.vtx - mesh.vtx)
norm(meshMsh.elt - mesh.elt)


disp('~~> Michto gypsilab !')



