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
%|    #    |   FILE       : nrtRayLabyrinthe.m                            |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray tracing with complex labyrinthe           |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nray = 1e5;
ord  = 7;
src  = [1.75 7 1];
mic  = [12 0 1];
rad  = 1;

% Initialize mesh
mesh = msh('labyrinthe.mesh');

% Initalize ray tracer
ray = ray(mesh,src,Nray);

% Plot mesh
figure
plot(mesh)
hold on
plot(ray)
axis equal;
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.3)
view(0,90)

% Measurement sphere
[X,Y,Z] = sphere(10);
surf(rad.*X+mic(1),rad.*Y+mic(2),rad.*Z+mic(3),ones(size(X)));
plot3(src(1),src(2),src(3),'pr','MarkerSize',20)

% Ray tracing
tic
ray = ray.tracer(ord);
toc

% Measures
tic
I = ray.measure(mic,rad);
toc

% Final ray
raym = ray.sub(I{end});
plot(raym)


disp('~~> Michto gypsilab !')


