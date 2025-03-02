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
%|    #    |   FILE       : nrtHmxCriticalDimension.m                     |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
%|  ( # )  |   CREATION   : 14.01.2019                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : Build fem H-Matrix of critical size           |
%|  `---'  |                                                              |
%+========================================================================+

% Nettoyage
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensions
Nx = 799 % (2^n)*100 - 1
Ny = 799

% Accuracy
tol = 1e-3;

% Meshes
Xmsh = mshSphere(Nx,1);
% Ymsh = mshCube(Ny,2*[1 1 1]);
Ymsh = Xmsh;%Ymsh.bnd;

% Domain
Xdom = dom(Xmsh,3);
Ydom = dom(Ymsh,3);

% Finite element
Xfem = fem(Xmsh,'P1');
Yfem = fem(Ymsh,'P1');

% Wave number or frequency (Hz)
k = 5;
f = (k*340)/(2*pi);

% Green kernel -> exp(1i*k*r)/r
green = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Particles charges
V = (-1+2*rand(length(Yfem),1)) + (-1+2i*rand(length(Yfem),1));

% Spatial representation of particles
figure
plot(Xmsh,'b')
hold on
plot(Ymsh,'r')
alpha(0.5)
axis equal 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCUL DIRECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FULL PRODUCT ~~~~~~~~~~~~~')
tic

% Full Matrix
tic
M = integral(Xdom,Ydom,Xfem,green,Yfem); 
toc

% Matrix-vector product
tic
ref = M*V;
toc

disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FAST PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FAST PRODUCT ~~~~~~~~~~~~~')

% Hierarchical Matrix
tic
Mh = integral(Xdom,Ydom,Xfem,green,Yfem,tol); 
toc

% Graphical representation
figure
spy(Mh)

% Matrix vector product
tic
sol = Mh * V;
toc

% Erreur
norm(ref-sol)/norm(ref)
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LU FACTORIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ LU FACTORIZATION ~~~~~~~~~~~~~')

% Factorization LU
tic
[Lh,Uh] = lu(Mh);
toc

% Graphical representation
figure
subplot(1,2,1); spy(Lh)
subplot(1,2,2); spy(Uh)

% Matrix vector product
tic
sol = Lh * (Uh * V);
toc

% Erreur
norm(ref-sol)/norm(ref)
disp(' ')


disp('~~> Michto gypsilab !')


