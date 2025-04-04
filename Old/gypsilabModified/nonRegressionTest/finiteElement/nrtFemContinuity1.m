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
%|    #    |   FILE       : nrtFemDirichletSquare.m                       |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Dirichlet condition with a square             |
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
Neig = 10;
L    = [1 0.5];
dl   = 1e-1;

% Original mesh
mesh  = mshSquare(Nvtx,L);
bound = mesh.bnd;
ctr   = mesh.ctr;

% Sub-mesh 1
mesh1        = mesh.sub(ctr(:,1)<=0);
mesh1.col(:) = 1;
int1         = setdiff(mesh1.bnd,bound);

% Sub mesh 2 with translation
mesh2          = mesh.sub(ctr(:,1)>=0);
mesh2.col(:)   = 2;
mesh2.vtx(:,1) = mesh2.vtx(:,1) + dl;
int2           = int1;
int2.vtx(:,1)  = int2.vtx(:,1) + dl;

% Final mesh
mesh  = union(mesh1,mesh2);
bound = mesh.bnd;

% Dirichlet
int = union(int1,int2);
dir = setdiff(bound,int);

% Graphical representation
figure
plot(mesh)
hold on
% plot(bound,'y')
plot(int1,'c')
plot(int2,'m')
plot(dir,'y')
axis equal

% Domain
omega = dom(mesh,7);

% Finites elements space
u = fem(mesh,'P2');

% Dirichlet condition
u = dirichlet(u,dir);

% Continuity
u = junction(u,int1,1,int2,-1);

% Mass matrix
tic
M = integral(omega,u,u);
toc

% Rigidity matrix
tic
K = integral(omega,grad(u),grad(u));
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
    axis equal off
    colorbar
end

% Analytical solutions of eigenvalues for an arbitrary cube
ref = zeros(Neig^2,1);
l = 1;
for i = 1:Neig
    for j = 1:Neig
        ref(l) = pi*sqrt( (i/L(1))^2 + (j/L(2))^2 );
        l = l+1;
    end
end
ref = sort(ref);
ref = ref(1:Neig);

% Error
sol = EV(1:Neig);
[ref sol abs(sol-ref)./ref]


disp('~~> Michto gypsilab !')


