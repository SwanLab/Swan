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
%|    #    |   FILE       : nrtFemDirichletString.m                       |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2018                                    |
%| ( === ) |   SYNOPSIS   : Dirichlet condition with a guitar string      |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nvtx = 100;
Neig = 10;
L    = 2;

% Meshes
mesh = mshSegment(Nvtx,L);

% Boundary
bound = mesh.bnd;

% Domain
omega = dom(mesh,3);
sigma = dom(bound,2);

% Finites elements space
u = fem(mesh,'P1');
v = fem(mesh,'P1');

% Dirichlet
u = dirichlet(u,bound);
v = dirichlet(v,bound);

% Graphical representation
plot(mesh,'w'); 
hold on
plot(bound,'r')
plot(omega)
plot(u,'go')
hold off
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.99)

% Mass matrix
tic
M = integral(omega,u,v);
toc
abs(sum(sum(M)) - 1)/1

% Rigidity matrix
tic
K = integral(omega,grad(u),grad(v));
toc
sum(sum(K))

% Find eigen values
tic
[V,EV] = eigs(K,M,2*Neig,'SM');
toc

% Normalization
V = V./(max(max(abs(V))));

% Sort by ascending order
[EV,ind] = sort(sqrt(real(diag(EV))));
[~,uni]  = unique(floor(1e2*EV));
V        = V(:,ind);

% Graphical representation
X = u.unk;
figure
for n = 1:9
    subplot(3,3,n)
    plot(X(:,1),V(:,n))
    title(['k = ',num2str(EV(n))])
    axis equal off
end

% Analytical solutions of eigenvalues for an arbitrary cube
ref = zeros(Neig,1);
for i = 1:Neig
    ref(i) = pi*sqrt( (i/L(1))^2 ) ;
end

% Error
sol = EV(1:Neig);
[ref sol abs(sol-ref)./ref]



disp('~~> Michto gypsilab !')



