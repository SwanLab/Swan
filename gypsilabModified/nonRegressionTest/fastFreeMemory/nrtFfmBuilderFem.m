%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2019.                             |
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
%|    #    |   FILE       : nrtFfmBuilderFem.m                            |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Non regression test for convolution using     |
%|  `---'  |                arbitrary kenel integral galerkin formulation |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ DEFINITIONS ~~~~~~~~~~~~~')

% Dimensions
Nx = 1e3
Ny = 2e3

% Accuracy
tol = 1e-3

% Meshes
Xmsh = mshSphere(Nx,1);
Ymsh = mshCube(Ny,2*[1 1 1]);
Ymsh = Ymsh.bnd;

% Domain
Xdom = dom(Xmsh,3);
Ydom = dom(Ymsh,3);

% Finite element
Xfem = fem(Xmsh,'P1');
Yfem = fem(Ymsh,'P1');

% Wave number or frequency (Hz)
stp = Xmsh.stp;
k   = 5

% Green kernel
green = '[exp(ikr)/r]'
Gxy   = @(X,Y) femGreenKernel(X,Y,green,k);

% Charges
V = (-1+2*rand(length(Yfem),1)) + (-1+2i*rand(length(Yfem),1));

% Spatial representation of particles
figure
plot(Xmsh,'b')
hold on
plot(Ymsh,'r')
alpha(0.5)
axis equal 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FULL PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FULL PRODUCT ~~~~~~~~~~~~~')

% Full Matrix
tic
M = integral(Xdom,Ydom,Xfem,Gxy,Yfem); 
toc

% Matrix-vector product
tic
ref = M * V;
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFM PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FFM PRODUCT ~~~~~~~~~~~~~')

% FFM
tic
Mv = integral(Xdom,Ydom,Xfem,green,k,Yfem,tol); 
toc

% FFM Matrix-vector product
tic
sol = Mv * V;
toc

% Error
norm(ref-sol)/norm(ref)



disp('~~> Michto gypsilab !')



