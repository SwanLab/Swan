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
%|    #    |   FILE       : nrtRaySphere.m                                |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray tracing with spherical mesh               |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nvtx = 1e4
Nray = 1e2
ord  = 5
rad  = 0.2;

% Spherical mesh
mesh = mshSphere(Nvtx,1);

% Plot mesh
figure
plot(mesh)
axis equal;
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.3)

% Initialize ray tracing
ray = ray(mesh,[0 0 0],Nray);

% Plot raytracing
plot(ray)

% Ray tracing
tic
ray = ray.tracer(ord);
toc

% Spherical measures
tic
[I,src] = ray.measure([0 0 0],rad);
toc

% Energy
src = cell2mat(src);
dst = sqrt(sum(src.^2,2));
dir = sparse(round(dst*10)+1,1,1);

% Error
nonzeros(dir)-Nray
figure
plot(full(dir))



disp('~~> Michto gypsilab !')



