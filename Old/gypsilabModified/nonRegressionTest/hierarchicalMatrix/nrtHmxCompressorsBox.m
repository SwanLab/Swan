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
%|    #    |   FILE       : nrtHmxCompressorBox.m                         |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Compare compressor for separated boxes        |
%|  `---'  |                                                              |
%+========================================================================+

clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Taille de la boite
edg = 1;

% Precision ou nombre de points de quadrature polynomiale
tol = 1e-3;

% Regular cube
n       = 10;
x       = 0:1/(n-1):1;
[x,y,z] = meshgrid(x,x,x);

% Recepteurs
Nx     = 1000;
% X      = edg*rand(Nx,3);
X      = edg*[x(:),y(:),z(:)];

% Emmeteurs
Ny     = 1000;
% Y      = edg*rand(Ny,3);
Y      = edg*[x(:),y(:),z(:)];
Y(:,1) = 2*edg + Y(:,1);

% Potentiel aleatoire aux emmeteurs
V  = (-1-1i) + (2+2i)*rand(Ny,1);

% Representation graphique
figure(1)
plot3(X(:,1),X(:,2),X(:,3),'*b',Y(:,1),Y(:,2),Y(:,3),'r*')
grid on
axis equal

% Noyau de green Helmholtz
k     = 5;
rxy   = @(X,Y) sqrt( (X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2 + (X(:,3)-Y(:,3)).^2 );
green = @(X,Y) exp(1i*k*rxy(X,Y))./rxy(X,Y);

% Gridding
[I,J] = ndgrid(1:Nx,1:Ny);

% Noyau de green
Gxy = green(X(I,:),Y(J,:));
Gxy = reshape(Gxy,Nx,Ny);
% Gxy = [Gxy zeros(Nx,Ny) ; ...
%     zeros(Nx,Ny) Gxy];
% Gxy = [Gxy zeros(Nx,2*Ny) ; ...
%     zeros(Nx,Ny) Gxy zeros(Nx,Ny) ; ...
%     zeros(Nx,2*Ny) Gxy];
% figure
% imagesc(Gxy)

% Calcul du rang
tic
rk = rank(Gxy,tol);
toc
rk

% Compression ACA, pivotage total
tic
[A,B] = hmxACA(Gxy,tol);
toc
size(A)
norm(A*B-Gxy,'fro')./norm(Gxy,'fro')

% Compression ACA, pivotage partiel
tic
[A,B] = hmxACA(X,Y,@(X,Y) green(X,Y),tol);
toc
size(A)
norm(A*B-Gxy,'fro')./norm(Gxy,'fro')

% % Recompression ACA, pivotage total
% tic
% [Ap,Bp] = hmxACA(A,B,tol);
% toc
% size(Ap)
% norm(Ap*Bp-Gxy,'fro')./norm(Gxy,'fro')

% Recompression RSVD
tic
[Ap,Bp,flag] = hmxRSVD(A,B,tol);
toc
size(Ap)
norm(Ap*Bp-Gxy,'fro')./norm(Gxy,'fro')

% Recompression QRSVD
tic
[Ap,Bp] = hmxQRSVD(A,B,tol);
toc
size(Ap)
norm(Ap*Bp-Gxy,'fro')./norm(Gxy,'fro')



disp('~~> Michto gypsilab !')


    


    




