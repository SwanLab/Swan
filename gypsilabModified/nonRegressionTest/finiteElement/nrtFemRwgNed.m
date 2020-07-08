%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Francois Alouges (c) 2017-2018.                            |
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
%|    #    |   FILE       : nrtFemRwgNed.m                                |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Francois Alouges                              |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Operators validation                          |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

disp('---- 2D problems --');
% Meshes
Nvtx = 1e3;
mesh = mshSphere(Nvtx,1);

% domains
omega = dom(mesh,3);

% FEM
u = fem(mesh,'RWG');
v = fem(mesh,'P1');
w = fem(mesh,'NED');

I = integral(omega, u, grad(v)) + integral(omega, div(u), v);
max(max(I))
J = integral(omega, w, nxgrad(v));
J = J - integral(omega, curl(w), v);
max(max(J))

disp('---- 3D problems --');
% Meshes
Nvtx = 1e3;
mesh = mshCube(Nvtx,[1 1 1]);

% domains
omega = dom(mesh,4);

% FEM
u = fem(mesh,'RWG');
v = fem(mesh,'P1');
w = fem(mesh,'NED');

f = {@(X) X(:,1),@(X) X(:,2),@(X) X(:,3)};
f = {@(X) ones(size(X,1),1),@(X) ones(size(X,1),1),@(X) ones(size(X,1),1)};
g = @(X) ones(size(X,1),1);
u0 = omega.interpolate(w,f);
I1 = integral(omega,f)
I2 = integral(omega,g,w);
I2{1}*u0
I2{2}*u0
I2{3}*u0
I = integral(omega,curl(w),curl(w))*u0;
max(abs(I))

clear all
% Meshes
Nvtx = 1e4;
mesh = mshCube(Nvtx,[1 1 1]);

% domains
omega = dom(mesh,4);
omega2 = dom(mesh,1);
% FEM
u = fem(mesh,'RWG');
w = fem(mesh,'NED');

disp('Masse RWG')
tic
I = integral(omega,u,u);
t1=toc
size(I)

disp('Rig RWG')
tic
I = integral(omega2,div(u),div(u));
t2=toc
size(I)

disp('Masse NED')
tic
I = integral(omega,w,w);
t3=toc
size(I)

disp('Rig NED')
tic
I = integral(omega2,curl(w),curl(w));
t4=toc
size(I)

disp('~~> Michto gypsilab !')


