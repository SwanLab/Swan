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
%|    #    |   FILE       : nrtHalfSquare.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Surfacic integration                          |
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
L    = [1 1];
tol  = 1e-3;

% Create square mesh
mesh = mshSquare(Nvtx,L);
plot(mesh)
hold on

% Colors half square
ctr = mesh.ctr;
mesh.col(ctr(:,1)>=0) = 1;

% Domain x>=0
omega    = dom(mesh.sub(mesh.col==1),3);
[Xqud,W] = omega.qud;

% Finite element space
u    = fem(mesh,'P2');
Xdof = u.dof;

% Graphical representation
plot(mesh)
hold on
plot(omega,'y.')
axis equal
alpha(0.99)
view(0,90)

% Dof to quadrature matrix
tic
M = u.uqm(omega);
toc

% /int_{[0,0.5],[-0.5,0.5]} 1 dx dy = 1 
(W' * M) * ones(size(Xdof,1),1) - 1/2

% /int_{[0,0.5],[-0.5,0.5]} x dx dy = 1/8 
(W' * M) * Xdof(:,1) - 1/8

% /int_{[0,0.5],[-0.5,0.5]} x dx dy = 1/8 
(W' * M) * (Xdof(:,1).^2) - 1/(3*8)

% Gradient
u  = grad(u);
Mp = u.uqm(omega);

% /int_{[0,0.5],[-0.5,0.5]} gradx(1) dx dy = 0 
(W' * Mp{1}) * ones(size(Xdof,1),1) - 0

% /int_{[0,0.5],[-0.5,0.5]} gradx(x) dx dy = 
(W' * Mp{1}) * Xdof(:,1) - 0.5

% /int_{[0,0.5],[-0.5,0.5]} gradx(x^2) dx dy = 1/8 
(W' * Mp{1}) * (Xdof(:,1).^2) - 1/4


% grad[psi]3
u   = grad(u,3);
tmp = u.uqm(omega);
norm(Mp{3}-tmp,'inf')

% N*[psi]
u   = ntimes(u);
tmp = u.uqm(omega);
norm(tmp{3}-M,'inf')

% N*[psi]3
u   = ntimes(u,3);
tmp = u.uqm(omega);
norm(tmp-M,'inf')

% Nxgrad[psi]
u   = nxgrad(u);
tmp = u.uqm(omega);
norm(tmp{2}-Mp{1},'inf')

% Nxgrad[psi]2
u   = nxgrad(u,2);
tmp = u.uqm(omega);
norm(tmp-Mp{1},'inf')


disp('~~> Michto gypsilab !')

