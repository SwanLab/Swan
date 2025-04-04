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
%|    #    |   FILE       : nrtFfmAlgebra.m                               |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2019                                    |
%| ( === ) |   SYNOPSIS   : Evaluate each ffm class function              |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Dimensions
Nx = 1e3;
Ny = 2e3;

% Accuracy
tol = 1e-3;

% Particles receptors X
mesh = mshSphere(Nx,1);
X    = mesh.vtx;

% Particles transmitters Y (=X or not)
mesh = mshCube(Ny,2*[1 1 1]);
Y    = mesh.vtx;
Ny   = size(Y,1);

% Type
type = 'double';

% Wave number or frequency (Hz)
k = 5;
f = (k*340)/(2*pi);

% Green kernel -> exp(1i*k*r)/r
green = '[exp(ikr)/r]';
Gxy   = @(X,Y) femGreenKernel(X,Y,green,k);

% Particles charges (multiples)
V = (-1+2*rand(Ny,2)) + (-1+2i*rand(Ny,2));

% Spatial representation of particles
figure
plot3(X(:,1),X(:,2),X(:,3),'*b',Y(:,1),Y(:,2),Y(:,3),'*r')
axis equal 


%%% Full matrix computation
disp('~~~~~~~~~~~~~ EXACT MATRIX AS REFERENCE~~~~~~~~~~~~~')
tic
M = zeros(Nx,Ny,type);
for i = 1:Nx
    M(i,:) = Gxy(X(i,:),Y).';
end
MV = M * V;
toc

tic
un = ones(Nx,1);
I  = spdiags([-un 2*un -un], [-Nx/2,0,Nx/2] , Nx, Ny);
IV = I*double(V);
toc
spy(I)
disp(' ')


%%% Fast & Furious computation
disp('~~~~~~~~~~~~~ FFM ~~~~~~~~~~~~~')
tic
Mv = ffm(X,Y,green,k,tol);
toc

tic
Mv2 = ffm(X,Y,green,k,0.1);
toc

disp(' ')


%%% Matrix-Vector product
disp('~~~~~~~~~~~~~ MATRIX * VECTOR ~~~~~~~~~~~~~')
tic
sol = Mv * V;
toc
ref = MV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = Mv2 * V;
toc
ref = MV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Scalar product
disp('~~~~~~~~~~~~~ SCALAR PRODUCT ~~~~~~~~~~~~~')
tic
tmp = (sqrt(2) .* Mv .* pi);
sol = tmp * V;
toc
ref = sqrt(2) * pi * MV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Uminus
disp('~~~~~~~~~~~~~ UMINUS ~~~~~~~~~~~~~')
tic
tmp = - Mv;
sol = tmp * V;
toc
ref = - MV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Addition
disp('~~~~~~~~~~~~~ ADDITION/SUBSTRACTION ~~~~~~~~~~~~~')
tic
tmp = 2.*M + Mv - M;
sol = tmp * V;
toc
ref = 2 * MV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = 2.*I + Mv - I;
sol = tmp * V;
toc
ref = MV + IV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Product
disp('~~~~~~~~~~~~~ MATRIX PRODUCT ~~~~~~~~~~~~~')
tic
tmp = Mv * (I' * Mv);
sol = tmp * V;
toc
ref = M * (I' * (M * V));
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = I' * Mv * (I'*I);
sol = tmp * V;
toc
ref = I' * (M * (I' * (I*V)));
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')



%%% Concatenation

%%% Single/Double


disp('~~> Michto gypsilab !')



