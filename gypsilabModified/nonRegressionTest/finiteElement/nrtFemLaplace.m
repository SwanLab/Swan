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
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtFemLaplace.m                               |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Solve Laplace equation for a disk             |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nvtx = 1e3;

% Meshes
mesh  = mshDisk(Nvtx,1);

% Domain
omega = dom(mesh,3);

% Finites elements space
u = fem(mesh,'P1');
v = fem(mesh,'P1');

% Graphical representation
plot(mesh,'w'); 
hold on
plot(omega)
plot(u,'go')
hold off
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.99)

% Rigidity matrix
tic
K = integral(omega,grad(u),grad(v));
toc
sum(sum(K))

% Mass matrix
tic
M = integral(omega,u,v);
toc
abs(sum(sum(M)) - pi)/pi

% Right hand side
f = @(X) X(:, 1).^2;
F = integral(omega,u,f);

% Solving
tic
uh = (K+M)\F ;
toc

% Plot the solution
figure
graph(u,uh);
title('Solution')
xlabel('X');   ylabel('Y');   zlabel('Z');
view(30,30)

disp('~~> Michto gypsilab !')


