%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Mathilde Boissier (c) 2017-2018.         |
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
%|    #    |   FILE       : elastcite2D.m                                 |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Mathilde Boissier           |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Elasticite lineaire 2D                        |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parametres
L      = [2 1];
N      = 1e4;
lambda = 2;
mu     = 3;

% Excitation
g = [50    0 0];
a = [1 -0.05 0];
b = [1 0.05  0];

% Maillage
mesh = mshSquare(N,L);

% Deplacement du point a
Rp             = sum( (mesh.vtx - ones(size(mesh.vtx,1),1)*a).^2 , 2);
[~,Ia]         = min(Rp);
mesh.vtx(Ia,:) = a;

% Deplacement du point b
Rp             = sum( (mesh.vtx - ones(size(mesh.vtx,1),1)*b).^2 , 2);
[~,Ib]         = min(Rp);
mesh.vtx(Ib,:) = b;

% Frontiere
bound = mesh.bnd;

% Maillage du bord gauche
ctr = bound.ctr;
I   = (ctr(:,1) == -1);
dir = bound.sub(I);

figure
plot(mesh)

% Maillage de l'excitation 
Ivtx = find ( (bound.vtx(:,1) >= a(1)) & (bound.vtx(:,2) >= a(2)) & ...
    (bound.vtx(:,1) <= b(1)) & (bound.vtx(:,2) <= b(2)) );
Ielt = [];
for i = 1:size(bound.elt,1)
   if sum(bound.elt(i,1) == Ivtx) && sum(bound.elt(i,2) == Ivtx)
        Ielt = [Ielt;i];
   end
end
neu = bound.sub(Ielt);

% Domaine d'integration
omega = dom(mesh,3);
sigma = dom(neu,2);

% Element finis
u   = fem(mesh,'P1');
phi = fem(mesh,'P1');

% Condition de Dirichlet
u   = dirichlet(u,dir);
phi = dirichlet(phi,dir);

% Affichage graphique
figure
plot(mesh)
hold on
plot(bound,'r')
plot(dir,'y')
plot(neu,'b')
grid on
axis equal

% Operateur e(phi):e(u)
A11 = integral(omega,grad(phi,1),grad(u,1)) + ...
    0.5 *  integral(omega,grad(phi,2),grad(u,2));

A12 = 0.5 *  integral(omega,grad(phi,2),grad(u,1));

A21 = 0.5 *  integral(omega,grad(phi,1),grad(u,2));

A22 = integral(omega,grad(phi,2),grad(u,2)) + ...
    0.5 *  integral(omega,grad(phi,1),grad(u,1));

% Operateur div(phi):div(u)
B11 = integral(omega,grad(phi,1),grad(u,1));
B12 = integral(omega,grad(phi,1),grad(u,2));
B21 = integral(omega,grad(phi,2),grad(u,1));
B22 = integral(omega,grad(phi,2),grad(u,2));

% Operateur vectoriel final
LHS = (2*mu) .* [A11 A12 ; A21 A22] + lambda .* [B11 B12 ; B21 B22];

% Representation graphique
figure
spy(LHS)

% Second membre
fx  = @(X) g(1) * ones(size(X,1),1);
fy  = @(X) g(2) * ones(size(X,1),1);
RHS = [integral(sigma,phi,fx) ; integral(sigma,phi,fy)];

% Solve probleme
X  = LHS \ RHS;
X1 = X(1:end/2);
X2 = X(end/2+1:end);

% Figure
figure
surf(phi,X1)
colorbar
axis equal

figure
surf(phi,X2)
colorbar
axis equal


disp('~~> Michto gypsilab !')



