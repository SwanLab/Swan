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
%|    #    |   FILE       : nrtFemDirichletCube.m                         |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Dirichlet condition with a cube               |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nvtx = 1e4;
Neig = 10;
L    = [1 0.5 0.5];

% Meshes
mesh  = mshCube(Nvtx,L);
% load tetmesh
% mesh = msh(X,tet);

% Boundary
bound = mesh.bnd;

% Domain
omega = dom(mesh,4);
sigma = dom(bound,3);

% Finites elements space
u = fem(mesh,'P1');
v = fem(mesh,'P1');

% Dirichlet
u = dirichlet(u,bound);
v = dirichlet(v,bound);

% Graphical representation
plot(mesh); 
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.1)
view(10,10)

% Mass matrix
tic
M = integral(omega,u,v);
toc

% Rigidity matrix
tic
K = integral(omega,grad(u),grad(v));
toc

% Find eigen values
tic
[V,EV] = eigs(K,M,2*Neig,'SM');
toc

% Normalization
V = V./(max(max(abs(V))));

% Sort by ascending order
[EV,ind] = sort(sqrt(real(diag(EV))));
V        = V(:,ind);

% Graphical representation
figure
for n = 1:9
    subplot(3,3,n)
    surf(u,V(:,n))
    title(['k = ',num2str(EV(n))])
    alpha(0.1)
    axis equal off
    view(10,10)
end

% Analytical solutions of eigenvalues for an arbitrary cube
ref = zeros(Neig^3,1);
l = 1;
for i = 1:Neig
    for j = 1:Neig
        for k = 1:Neig
            ref(l) = pi*sqrt( (i/L(1))^2 + (j/L(2))^2 +  (k/L(3))^2 );
            l = l+1;
        end
    end
end
ref = sort(ref);
ref = ref(1:Neig);

% Error
sol = EV(1:Neig);
[ref sol abs(sol-ref)./ref]



disp('~~> Michto gypsilab !')


