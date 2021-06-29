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
%|    #    |   FILE       : nrtHelmholtzSDrad.m                           |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
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
N   = 1e3
tol = 1e-3
typ = 'P1'
gss = 3
X0  = 0*[0.1 0.1 0.1]
k   = 5

% Boundary mesh
mesh = mshSphere(N,1);

% Radiative mesh
radiat     = mshSquare(10*N,[5 5]);
radiat.vtx = [radiat.vtx(:,1) zeros(size(radiat.vtx,1),1) radiat.vtx(:,2)];

% Domain
sigma = dom(mesh,gss);    

% Finite element
u = fem(mesh,typ);

% Mesh representation
figure
plot(mesh)
hold on
plot(radiat)
hold off
axis equal
axis(2.5*[-1 1 -1 1 -1 1])
view(0,0)

% Frequency adjusted to maximum esge size
stp = mesh.stp;
if (k > 1/stp(2))
    warning('frequency is too hight for mesh')
end

% Spherical sources
SW      = @(X) femGreenKernel(X,X0,'[exp(ikr)/r]',k);
dxSW{1} = @(X) femGreenKernel(X,X0,'gradx[exp(ikr)/r]1',k);
dxSW{2} = @(X) femGreenKernel(X,X0,'gradx[exp(ikr)/r]2',k);
dxSW{3} = @(X) femGreenKernel(X,X0,'gradx[exp(ikr)/r]3',k);

% Jump trace : mu = [p] = pi - pe    lambda = [dnp] = dnpi - dnpe 
mu     = - integral(sigma,u,SW);    % pi = 0
lambda = - integral(sigma,ntimes(u),dxSW);  % dnpi = 0

% Mass matrix
Id = integral(sigma,u,u);

% Single Layer
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
S   = 1/(4*pi) .* ( integral(radiat.vtx,sigma,Gxy,u)  ...
    + regularize(radiat.vtx,sigma,'[1/r]',u) );

% Double layer
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);
D        = 1/(4*pi) .* ( integral(radiat.vtx,sigma,dyGxy,ntimes(u)) ...
    + regularize(radiat.vtx,sigma,'grady[1/r]',ntimes(u)) );

% Radiation with integral representation
sol = S * (Id\lambda) - D * (Id\mu);
% sol = S * (Id\lambda);
% sol = - D * (Id\mu);

% Analytic solution : Gk(X,X0) exterieur et 0 interieur
ref = SW(radiat.vtx);
ind = find(sqrt(sum(radiat.vtx.^2,2)) > 1);
norm(ref(ind)-sol(ind))/norm(ref(ind))

% Graphical representation
theta = 0:2*pi/1e3:2*pi;
figure
plot(radiat,real(ref))
hold on
plot3(cos(theta),zeros(size(theta)),sin(theta),'k--','LineWidth',3)
hold off
title('Analytic solution')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
colorbar
caxis([-1 1])
axis equal
axis(2.5*[-1 1 -1 1 -1 1])
view(0,0)

figure
plot(radiat,real(sol))
hold on
plot3(cos(theta),zeros(size(theta)),sin(theta),'k','LineWidth',3)
hold off
title('Radiation solution')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
colorbar
caxis([-1 1])
axis equal
axis(2.5*[-1 1 -1 1 -1 1])
view(0,0)



disp('~~> Michto gypsilab !')


