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
%|    #    |   FILE       : nrtBmmAlgebra.m                               |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Block matrix algebra                          |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
run('../../addpathGypsilab')

% Dimensions
N = 200;

% Accuracy
tol = 1e-3;

% Particles receptors X
mesh  = mshSphere(N,1);
omega = dom(mesh,3);
phi0  = fem(mesh,'P0');
phi1  = fem(mesh,'P1');

% Wave number or frequency (Hz)
k = 5;

% Green kernel -> exp(1i*k*r)/r
green = '[exp(ikr)/r]';
Gxy   = @(X,Y) femGreenKernel(X,Y,green,k) + ...
    sqrt(N).*(X(:,1)==Y(:,1)).*(X(:,2)==Y(:,2)).*(X(:,3)==Y(:,3))  ;

% Particles charges (multiples)
V0 = (-1+2*rand(length(phi0),2)) + (-1+2i*rand(length(phi0),2));
V1 = (-1+2*rand(length(phi1),2)) + (-1+2i*rand(length(phi1),2));

% Spatial representation of particles
figure
plot(mesh)
axis equal 

% All forms
tic
Ah = integral(omega,omega,phi0,Gxy,phi0,tol);
Bv = integral(omega,omega,phi0,green,k,phi1,tol);
Cs = integral(omega,phi1,phi0);
Df = integral(omega,omega,phi1,Gxy,phi1);
toc

%%% Single Builder
disp('~~~~~~~~~~~~~ SINGLE BUILDERS ~~~~~~~~~~~~~')
Mb  = bmm(Ah);
sol = Mb * V0;
ref = Ah * V0;
norm(ref-sol,'inf')/norm(ref,'inf')

Mb  = bmm(Bv);
sol = Mb * V1;
ref = Bv * V1;
norm(ref-sol,'inf')/norm(ref,'inf')

Mb  = bmm(Cs);
sol = Mb * V0;
ref = Cs * V0;
norm(ref-sol,'inf')/norm(ref,'inf')

Mb  = bmm(Df);
sol = Mb * V1;
ref = Df * V1;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Multi Builder
disp('~~~~~~~~~~~~~ MULTI BUILDERS ~~~~~~~~~~~~~')
Mb  = bmm({Ah,Cs';Cs,Df});
Mr  = [full(Ah) Cs' ; Cs Df];
Vb  = [V0;V1];
sol = Mb * Vb;
ref = Mr * Vb;
norm(ref-sol,'inf')/norm(ref,'inf')

figure
spy(Mb)

disp(' ')


%%% Full conversion
disp('~~~~~~~~~~~~~ FULL CONVERSION ~~~~~~~~~~~~~')
sol = full(Mb) * Vb;
ref = Mr * Vb;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Sparse conversion
disp('~~~~~~~~~~~~~ SPARSE CONVERSION ~~~~~~~~~~~~~')
sol = sparse(Mb) * Vb;
ref = Mr * Vb;
norm(ref-sol,'inf')/norm(ref,'inf')

figure
spy(sparse(Mb))

disp(' ')


%%% Transposition
disp('~~~~~~~~~~~~~ TRANSPOSITION ~~~~~~~~~~~~~')
sol = Vb.' * Mb.';
ref = (Mr * Vb).';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Transposition
disp('~~~~~~~~~~~~~ CONJUGATE TRANSPOSITION ~~~~~~~~~~~~~')
sol = Vb' * Mb';
ref = (Mr * Vb)';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Concatenation
disp('~~~~~~~~~~~~~ CONCATENATION ~~~~~~~~~~~~~')
Ab  = Mb;
Bb  = bmm({Bv;Df});
Cb  = bmm({Ah,Bv});
Db  = bmm(Cs');
Mb2 = [Ab , Bb ; Cb , Db];
sol = Mb2 * [Vb;V1];
ref = [Ab*Vb + Bb*V1;Cb*Vb+Db*V1];
norm(ref-sol,'inf')/norm(ref,'inf')

spy(Mb2)

disp(' ')


%%% Scalar product
disp('~~~~~~~~~~~~~ SCALAR PRODUCT ~~~~~~~~~~~~~')
sol = (-(3.*Mb.*2i)) * Vb;
ref = (-3*2i) * (Mr*Vb);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Addition
disp('~~~~~~~~~~~~~ ADDITION ~~~~~~~~~~~~~')
sol = (Mb+Mb) * Vb;
ref = 2*(Mr*Vb);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Multiplication
disp('~~~~~~~~~~~~~ MULTIPLICATION ~~~~~~~~~~~~~')
sol = (Mb.'*Mb) * Vb;
ref = Mr.'*(Mr*Vb);
norm(ref-sol,'inf')/norm(ref,'inf')

figure
spy(Mb.'*Mb)

disp(' ')


%%% LU factorization
disp('~~~~~~~~~~~~~ LU FACTORISATION ~~~~~~~~~~~~~')
[Lb,Ub] = lu(Mb);
sol = Lb * (Ub * Vb);
ref = Mb * Vb;
norm(ref-sol,'inf')/norm(ref,'inf')

figure
subplot(1,2,1); spy(Lb)
subplot(1,2,2); spy(Ub)

disp(' ')


%%% Solve LU
disp('~~~~~~~~~~~~~ SOLVE LU SYSTEM ~~~~~~~~~~~~~')
sol = Mb \ Vb;
ref = Mr \ Vb;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Inversion
disp('~~~~~~~~~~~~~ INVERSION ~~~~~~~~~~~~~')
Mbm1 = inv(Mb);
sol  = Mbm1 * Vb;
ref  = Mr \ Vb;
norm(ref-sol,'inf')/norm(ref,'inf')

figure
spy(Mbm1)

disp(' ')




disp('~~> Michto gypsilab !')



