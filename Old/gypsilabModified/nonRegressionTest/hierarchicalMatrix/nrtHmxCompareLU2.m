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
%|    #    |   FILE       : nrtHmxCompareLU.m                             |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Compar sparse LU and H-Matrix sparse LU       |
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
L    = [1 0.1 0.1];
tol  = 1e-6;

% Cube mesh
mesh = mshCube(Nvtx,L);

% Reorder
% mesh = symrcm(mesh);

% Domain
omega = dom(mesh,4);

% Graphical representation
plot(mesh);
hold on
plot(omega); 
axis equal;
title('Domain representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.1)
hold on

% Finites elements space
u = fem(mesh,'P1');
v = fem(mesh,'P1');

% Random vector at dof
V = rand(size(v.dof,1),1);

% Mass matrix
tic
M = integral(omega,u,v);
toc
ref = M * V;
figure
spy(M)
drawnow

% LU factorisation
tic
[L,U] = lu(M);
toc
norm(L*(U*V) - ref)/norm(ref)
figure
spy(L)

% Cholesky
tic
Uc = chol(M);
toc

% Error
ref = U\(L\V);
sol = Uc\(Uc'\V);
norm(ref-sol,'inf')/norm(ref,'inf')

% H-Matrix inverse
Mh = hmx(u.unk,u.unk,M,tol);
tic
Mhm1 = inv(Mh);
toc

% Error
ref = M\V;
sol = Mhm1*V;
norm(ref-sol,'inf')/norm(ref,'inf')

figure
subplot(1,2,1)
spy(Mh)
subplot(1,2,2)
spy(Mhm1)

% Inverse using facto
N     = size(M,1)
I     = symrcm(M)';
[L,U] = lu(M(I,I));
ind   = [(1:N)',zeros(N,2)];
Lh    = hmx(ind,ind,L,tol)
Uh    = hmx(ind,ind,U,tol)
Ih    = hmx(ind,ind,speye(size(M)),tol)
tic
Mhm1 = Uh\(Lh\Ih) 
toc

close all
figure
spy(Mhm1)




% Error
ref = M\V;
sol(I) = U\(L\V(I));
norm(ref-sol,'inf')/norm(ref,'inf')
sol(I) = Uh\(Lh\V(I));
norm(ref-sol,'inf')/norm(ref,'inf')
sol(I) = Mhm1*V(I);
norm(ref-sol,'inf')/norm(ref,'inf')




pouetPouet

% % Cuthill mckee
% I = symrcm(M);
% figure
% spy(M(I,I))

% H-Matrix
tic
Mh = hmx(u.dof,v.dof,M,1e-3);
toc
norm(Mh*V - ref)/norm(ref)
figure
spy(Mh)

% LhUh factorisation
tic
[Lh,Uh] = lu(Mh);
toc
norm(Lh*(Uh*V) - ref)/norm(ref)
figure
spy(Lh)
drawnow

% Solver
tic
Mm1V = U \ (L \ V);
toc
tic
Mhm1V = Uh \ (Lh \ V);
toc
tic
ref = M \ V;
toc
norm(Mm1V - ref)/norm(ref)
norm(Mhm1V - ref)/norm(ref)

whos



disp('~~> Michto gypsilab !')