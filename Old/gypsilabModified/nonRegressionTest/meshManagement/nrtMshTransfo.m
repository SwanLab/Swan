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
%|    #    |   FILE       : nrtMshTransfo.m                               |
%|    #    |   VERSION    : 0.52                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF : 21.06.2019                                    |
%| ( === ) |   SYNOPSIS   : Elementar geometry operations                 |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Create mesh
Nvtx = 1e3;
L    = [1 1];
mesh = mshSquare(Nvtx,L);

% Translation
mesh2 = translate(mesh,[1 1 1]);
mesh2.col(:) = 1;

% Rotation
mesh3 = mesh;
mesh3 = rotate(mesh3,[1 0 0],pi/4);
mesh3 = rotate(mesh3,[0 1 0],pi);
mesh3 = rotate(mesh3,[0 0 1],pi/2);
mesh3.col(:) = 2;

% Graphical representation
figure
plot(mesh)
axis equal
hold on
plotNrm(mesh,'r')
plot(mesh2)
plot(mesh3)
plotNrm(mesh3,'r')
view(45,45)
xlabel('X')
ylabel('Y')
zlabel('Z')

% Planar split
mesh          = mshSphere(Nvtx,1);
[mesh1,mesh2] = split(mesh,[0 0 0.7],[1 1 1]);

% Graphical representation
figure
hold on
plot(mesh1,'r')
plot(mesh2,'y')
axis equal
view(45,45)
xlabel('X')
ylabel('Y')
zlabel('Z')




disp('~~> Michto gypsilab !')
