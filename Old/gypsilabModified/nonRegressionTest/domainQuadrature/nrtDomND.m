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
%|    #    |   FILE       : nrtDomND.m                                    |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Multi-dimensional quadrature                  |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Tetra mesh
Nvtx = 1e2;
L    = [5 4 3];
mesh = mshCube(Nvtx,L);

% Boundary mesh
bound = mesh.bnd;
ctr   = bound.ctr;
bound.col(ctr(:,3)==0.5*L(3)) = 1;

% Upperface
upperface = bound.sub(bound.col==1); 
rectangle = upperface.bnd;

% Domain
omega = dom(mesh,5);
figure
plot(mesh)
hold on
plot(omega)
alpha(0.2)
axis equal
view(30,30)

% Boundary domain
sigma = dom(bound,3);
figure
plot(bound)
hold on
plot(sigma,'.y')
plotNrm(sigma,'y')
alpha(0.5)
axis equal
view(30,30)

% Rectangle
gamma = dom(rectangle,3);
figure
plot(rectangle)
hold on
plot(gamma)
alpha(0.5)
axis equal
view(30,30)

% Numerical integration on the cube
fct = @(X) ones(size(X,1),1);
V   = integral(omega,fct);
norm(V-prod(L),'inf')

% Numerical integration on the bundary of the cube
fct = @(X) ones(size(X,1),1);
S   = integral(sigma,fct);
norm(S-2*(L(1)*L(2)+L(2)*L(3)+L(1)*L(3)),'inf')

% Numerical integration on the rectangle
fct = @(X) ones(size(X,1),1);
P   = integral(gamma,fct);
norm(P-2*(L(1)+L(2)),'inf')


disp('~~> Michto gypsilab !')


