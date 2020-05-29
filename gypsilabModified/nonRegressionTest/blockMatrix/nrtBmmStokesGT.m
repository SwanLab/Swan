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
%|    #    |   FILE       : nrtBmmStokesGT.m                              |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 31.10.2018                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
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
N   = 256
x0  = [1 2 0.5];
rad = [5 3 2];
err = zeros(2,length(N));
h   = zeros(size(N));
tol = 1e-6;

% Mesh unit sphere
mesh = mshSphere(N,1);

% Sphere to ellipsoid
mesh.vtx = mesh.vtx .* (ones(N,1)*rad);

% Quadrature
gamma = dom(mesh,3);

% Finite element
phi = fem(mesh,'P1');
unk = phi.unk;
un  = ones(length(phi),1);

% Mass matrix
I = integral(gamma,phi,phi);
Z = zeros(size(I));
I = [I Z Z ; Z I Z ; Z Z I];

% Use of operator symetry
ind  = [1 1 ; 1 2 ; 1 3 ; 2 2 ; 2 3 ; 3 3];
ind  = [kron(ones(4,1),ind) kron((0:3)',ones(6,1))];
M    = cell(24,1);
C    = cell(24,1);

% Loop for operators
tic
parfor l = 1:length(ind)
    % Local indices
    i = ind(l,1);
    j = ind(l,2);
    k = ind(l,3);
    
    % Stokeslet
    if (k == 0)
        name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
        green = @(X,Y) femGreenKernel(X,Y,name,[]);
        M{l}  = 1/(8*pi) .* (integral(gamma,gamma,phi,green,phi,tol) + ...
            regularize(gamma,gamma,phi,name,phi));
    
    else
        % Stresslet
        name  = ['[rirjrk/r^5]',num2str(i),num2str(j),num2str(k)];
        green = @(X,Y) femGreenKernel(X,Y,name,[]);
        M{l}  = -6/(8*pi) .* integral(gamma,gamma,phi,green,ntimes(phi,k),tol);
        
        % Correction
        P    = -6/(8*pi) .* integral(gamma.qud,gamma,green,ntimes(phi,k),tol) * un;
        C{l} = integral(gamma,phi,@(X)P,phi);
    end
end
toc

% Block-Matrix operators
G = bmm({M{1}    M{2}    M{3}    ; M{2}    M{4}    M{5}    ; M{3}    M{5}    M{6}});
T = bmm({M{6+1}  M{6+2}  M{6+3}  ; M{6+2}  M{6+4}  M{6+5}  ; M{6+3}  M{6+5}  M{6+6}}) + ...
    bmm({M{12+1} M{12+2} M{12+3} ; M{12+2} M{12+4} M{12+5} ; M{12+3} M{12+5} M{12+6}}) + ...
    bmm({M{18+1} M{18+2} M{18+3} ; M{18+2} M{18+4} M{18+5} ; M{18+3} M{18+5} M{18+6}});

% Block-Matrix correction
C = bmm({C{6+1}  C{6+2}  C{6+3}  ; C{6+2}  C{6+4}  C{6+5}  ; C{6+3}  C{6+5}  C{6+6}}) + ...
    bmm({C{12+1} C{12+2} C{12+3} ; C{12+2} C{12+4} C{12+5} ; C{12+3} C{12+5} C{12+6}}) + ...
    bmm({C{18+1} C{18+2} C{18+3} ; C{18+2} C{18+4} C{18+5} ; C{18+3} C{18+5} C{18+6}});

% Graphical representation
figure
spy(G)
figure
spy(T)
figure
spy(C)

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
mu     = cell2mat(mu);
lambda = cell2mat(lambda);
toc

% Analytic solution : Gi1
ref = - mu;

% Stokes radiation without correction
sol = -G*(I\lambda) + T*(I\mu) - 0.5*mu;
norm(ref-sol)/norm(ref)

% Stokes radiation with correction
sol = -G*(I\lambda) + T*(I\mu) - C*(I\mu);
norm(ref-sol)/norm(ref)



disp('~~> Michto gypsilab !')



