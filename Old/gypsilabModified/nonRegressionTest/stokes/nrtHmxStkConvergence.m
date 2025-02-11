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
%|    #    |   FILE       : nrtHmxStkConvergence.m                        |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 31.10.2018                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Stoke Galerkin of an ovoid using H-Matrix with|
%|  `---'  |                stokeslet G and stresslet T                   |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N   = 2.^(6:8)
x0  = [1 2 0.5];
rad = [5 3 2];
err = zeros(2,length(N));
h   = zeros(size(N));
tol = 1e-6;

% Loop over mesh
for n = 1:length(N)
    % Mesh unit sphere
    mesh = mshSphere(N(n),1);
    
    % Sphere to ellipsoid
    mesh.vtx = mesh.vtx .* (ones(N(n),1)*rad);
    
    % Mesh step (mean)
    tmp  = mesh.stp;
    h(n) = tmp(3);
    
    % Quadrature
    gamma = dom(mesh,3);
    
    % Finite element
    phi = fem(mesh,'P1');
    unk = phi.unk;
    un  = ones(length(phi),1);

    % Mass matrix
    M = integral(gamma,phi,phi);
    
    % Use of operator symetry
    ind = [1 1 ; 2 2 ; 3 3 ; 1 2 ; 1 3 ; 2 3];
    G   = cell(6,1);
    T   = cell(6,1);
    C   = cell(6,1);
    I   = cell(6,1);
    
    % Loop for operators
    tic
    for l = 1:length(ind)
        % Local indices
        i = ind(l,1);
        j = ind(l,2);
        
        % Using r = x - y
        % Single layer : G = \int_gamma \int_gamma  1/(8pi) (\delta_ij/r + r_i*r_j/|r|^3)
        name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
        green = @(X,Y) 1/(8*pi) .* femGreenKernel(X,Y,name,[]);
        G{l}  = integral(gamma,gamma,phi,green,phi,tol);
        
        % Regularization
        G{l} = G{l} + 1/(8*pi) .* regularize(gamma,gamma,phi,name,phi);
        
        % Double layer  : T = \int_gamma \int_gamma -6/(8pi) (r_i*r_j*(r.n)/|r|^5)
        % Double lumped : P = \int_gamma -6/(8pi) (r_i*r_j*(r.n)/|r|^5)
        T{l} = zeros(G{l});
        P    = 0;
        for k = 1:3
            name  = ['[rirjrk/r^5]',num2str(i),num2str(j),num2str(k)];
            green = @(X,Y) -6/(8*pi) .* femGreenKernel(X,Y,name,[]);
            T{l}  = T{l} + integral(gamma,gamma,phi,green,ntimes(phi,k),tol);
            P     = P + integral(gamma.qud,gamma,green,ntimes(phi,k),tol) * un;
        end
        
        % Correction
        C{l} = integral(gamma,phi,@(X)P,phi);
        
        % Identity matrix
        if (i == j)
            I{l} = M;
        else
            I{l} = sparse(size(M,1),size(M,2));
        end
    end
    toc
    
    % Right hand side
    tic
    mu     = cell(3,1);
    lambda = cell(3,1);
    for i = 1:3
        % mu = [u] = u_int - u_ext = - Gi1 = - Gi1(x0,y)
        name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(1)];
        green = @(Y) 1/(8*pi) .* femGreenKernel(x0,Y,name,[]);
        mu{i} = - integral(gamma,phi,green);
        
        % lambda = [sigma] = - Ti1
        lambda{i} = 0;
        for k = 1:3
            name      = ['[rirjrk/r^5]',num2str(i),num2str(1),num2str(k)];
            green     = @(Y) -6/(8*pi) .* femGreenKernel(x0,Y,name,[]);
            lambda{i} = lambda{i} - integral(gamma,ntimes(phi,k),green);
        end
    end
    toc
    
    % Convert cells to H-Matrix and vectors
    tic
    G      = [G{1} G{4} G{5} ; G{4} G{2} G{6} ; G{5} G{6} G{3}];
    T      = [T{1} T{4} T{5} ; T{4} T{2} T{6} ; T{5} T{6} T{3}];
    C      = [C{1} C{4} C{5} ; C{4} C{2} C{6} ; C{5} C{6} C{3}];
    I      = [I{1} I{4} I{5} ; I{4} I{2} I{6} ; I{5} I{6} I{3}];
    mu     = cell2mat(mu);
    lambda = cell2mat(lambda);
    toc
    
%     % Validation comparing to BEM
%     norm(full(G))
%     norm(full(T))
%     norm(full(C))
%     norm(full(I))
%     norm(mu)
%     norm(lambda)
    
    % Analytic solution : Gi1
    ref = - mu;
        
    % Stokes radiation without correction
    sol = -G*(I\lambda) + T*(I\mu) - 0.5*mu;

    % Error L2
    tmp      = ref-sol;
    err(1,n) = sqrt(tmp'*(I\tmp));
    norm(tmp)./norm(ref)

    % Stokes radiation with correction
    sol = -G*(I\lambda) + T*(I\mu) - C*(I\mu);

    % Relative error L2 and inf
    tmp      = ref-sol;
    err(2,n) = sqrt(tmp'*(I\tmp));
    norm(tmp)./norm(ref)
    '========================================'
end

% Graphical representation
figure
loglog(h,h*3e-3 ,'--b')
hold on
loglog(h,h.^2*1e-3 ,'--r')
loglog(h,err(1,:),'-+b')
loglog(h,err(2,:),'-+r')
grid on
xlabel('h')
ylabel('Error L2')
legend({'Slope 1','Slope 2','T uncorr','T corr'})
title('Stokes - Fibonacci ellipsoide')
hold off




disp('~~> Michto gypsilab !')



