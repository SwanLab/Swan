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
%|    #    |   FILE       : nrtStkRadiation.m                             |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 31.10.2018                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Stoke radiation of a unit sphere using BEM    |
%|  `---'  |                with stokeslet G and stresslet T              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N  = 300;
n  = 100;
x0 = [0.1 0.2 0.3];

% Mesh unit sphere
mesh = mshSphere(N,1);

% Quadrature
gamma = dom(mesh,3);

% Finite element and unknowns
phi = fem(mesh,'P1');
unk = phi.unk;

% Radiation particles
mesh2 = mshSphere(n,2);
X     = mesh2.vtx;

% Graphical rep
figure
plot(mesh)
hold on
plotNrm(mesh)
plot(gamma)
plot(phi,'*r')
plot(msh(X))
alpha(0.5)
axis equal

% Initalization of block matrix for each coordinate (x,y,z)
G      = cell(3,3);
T      = cell(3,3);
mu     = cell(3,1);
lambda = cell(3,1);

% Loop for each coordinate
for i = 1:3
    for j = 1:3
        % Using r = x - y
        % Single layer : G = \int_gamma 1/(8pi) (\delta_ij/r + r_i*r_j/|r|^3) 
        name   = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
        green  = @(X,Y) 1/(8*pi) .* femGreenKernel(X,Y,name,[]);
        G{i,j} = integral(X,gamma,green,phi);
        
        % Double layer : T = \int_gamma -6/(8pi) (r_i*r_j*(r.n)/|r|^5)
        T{i,j} = 0;
        for k = 1:3
            name   = ['[rirjrk/r^5]',num2str(i),num2str(j),num2str(k)];
            green  = @(X,Y) -6/(8*pi) .* femGreenKernel(X,Y,name,[]);
            T{i,j} = T{i,j} + integral(X,gamma,green,ntimes(phi,k));
        end
    end
    
    % mu = [u] = u_int - u_ext = - Gi1 = - Gi1(x0,y)
    name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(1)];
    mu{i} = -1/(8*pi) * femGreenKernel(x0,phi.unk,name,[]);
    
    % lambda = [sigma] = - Ti1, with nk = xk (because boundary is unit sphere)
    lambda{i} = 0;
    for k = 1:3
        name      = ['[rirjrk/r^5]',num2str(i),num2str(1),num2str(k)];
        lambda{i} = lambda{i} - (-6/(8*pi)) * femGreenKernel(x0,phi.unk,name,[]) .* unk(:,k);
    end
end

% Convert cells to full matrix
G      = cell2mat(G);
T      = cell2mat(T);
mu     = cell2mat(mu);
lambda = cell2mat(lambda);

% Stokes radiation in domain : ui(x) = - sum_j \int_gamma Gij(x,y) lambda_j dy ...
%              + sum_j \int_gamma Tij(x,y).n(y) mu_j dy 
sol = -G*lambda + T*mu ; 

% Analytic solution : Gi1
ref = cell(3,1);
for i = 1:3
    name   = ['[ij/r+rirj/r^3]',num2str(i),num2str(1)];
    ref{i} = 1/(8*pi) * femGreenKernel(X,x0,name,[]);
end
ref = cell2mat(ref);

% Relative error L2 and inf
norm(ref-sol)/norm(ref)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation of solution
figure(2)
for i = 1:3
    subplot(1,3,i)
    ind = (i-1)*n + (1:n);
    plot(mesh2,sol(ind))
    axis equal
    colorbar
    title(['Component ',num2str(i)])
end





disp('~~> Michto gypsilab !')




