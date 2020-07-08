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
%|    #    |   FILE       : nrtHelmholtz2dSDrad.m                         |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Martin Averseng             |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : Solve dirichlet scatering problem with single |
%|  `---'  |                layer potential                               |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N   = 5e2
typ = 'P1'
gss = 3
X0  = zeros(1,3)
k   = 5

% Boundary mesh
mesh = mshCircle(N,1);

% Ellipse
% mesh.vtx(:,1) = 0.7*mesh.vtx(:,1);
% mesh.vtx(:,2) = 0.5*mesh.vtx(:,2);

% Radiative mesh
radiat = mshSquare(10*N,[5 5]);

% Domain
sigma = dom(mesh,gss);

% Finite element
u = fem(mesh,typ);

% Mesh representation
figure
plot(mesh)
hold on
plotNrm(mesh)
plot(radiat,'w')
plot(sigma)
plot(u)
hold off
axis equal
axis(2.5*[-1 1 -1 1 -1 1])
alpha(0.5)
view(0,90)

% Spherical sources
SW      = @(X) femGreenKernel(X,X0,'[H0(kr)]',k);
dxSW{1} = @(X) femGreenKernel(X,X0,'gradx[H0(kr)]1',k);
dxSW{2} = @(X) femGreenKernel(X,X0,'gradx[H0(kr)]2',k);
dxSW{3} = @(X) femGreenKernel(X,X0,'gradx[H0(kr)]3',k);

% Green kernels
Gxy      = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k);
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]1',k);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]2',k);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]3',k);

% Jump trace : mu = [p] = pi - pe    lambda = [dnp] = dnpi - dnpe 
mu     = - integral(sigma,u,SW);    % pi = 0
lambda = - integral(sigma,ntimes(u),dxSW);  % dnpi = 0

% Validation for mu
sol = -sum(mu);
ref = besselh(0,k) * 2*pi;
norm(ref-sol)/norm(ref)

% Validation for lambda
sol = -sum(lambda);
ref = - k*besselh(1,k) * 2*pi;
norm(ref-sol)/norm(ref)

% Mass matrix
Id = integral(sigma,u,u);

% Single Layer
tic
S   = 1i/4 .* integral(radiat.vtx,sigma,Gxy,u);
Sr  = -1/(2*pi) .* regularize(radiat.vtx,sigma,'[log(r)]',u);
S   = S + Sr;
toc

% Double layer
tic
D  = 1i/4 .* integral(radiat.vtx,sigma,dyGxy,ntimes(u));
Dr = -1/(2*pi) .* regularize(radiat.vtx,sigma,'grady[log(r)]',ntimes(u));
D  = D + Dr;
toc

% Radiation with integral representation
sol = S * (Id\lambda) - D * (Id\mu);

% Analytic solution : Gk(X,X0) exterieur et 0 interieur
ref = SW(radiat.vtx);
ind = find(sqrt(sum(radiat.vtx.^2,2)) >= 1.01);
norm(ref(ind)-sol(ind))/norm(ref(ind))
norm(ref(ind)-sol(ind),'inf')/norm(ref(ind),'inf')

% Graphical representation
theta = 0:2*pi/1e3:2*pi;
figure
plot(radiat,real(ref))
hold on
plot(mesh,'k')
hold off
title('Analytic solution')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
colorbar
caxis([-1 1])
axis equal
axis(2.5*[-1 1 -1 1 -1 1])
view(0,90)

figure
plot(radiat,real(sol))
hold on
plot(mesh,'k')
hold off
title('Radiation solution')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
colorbar
caxis([-1 1])
axis equal
axis(2.5*[-1 1 -1 1 -1 1])
view(0,90)

% Single Layer
tic
S  = 1i/4 .* integral(sigma,sigma,u,Gxy,u);
Sr = -1/(2*pi) .* regularize(sigma,sigma,u,'[log(r)]',u);
S  = S + Sr;
toc

% Double layer
tic
D  = 1i/4 .* integral(sigma,sigma,u,dyGxy,ntimes(u));
Dr = -1/(2*pi) .* regularize(sigma,sigma,u,'grady[log(r)]',ntimes(u));
D  = D + Dr;
toc

% Radiation with integral representation
sol = S * (Id\lambda) - D * (Id\mu) - 0.5*mu;

% Analytic solution : Gk(X,X0) exterieur et 0 interieur
ref = - mu;
norm(ref-sol)/norm(ref)
norm(ref-sol,'inf')/norm(ref,'inf')




disp('~~> Michto gypsilab !')



