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
%|    #    |   FILE       : nrtHmxAlgebra.m                               |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Evaluate each hmx class function              |
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
Ny = 1e3;

% Accuracy
tol = 1e-3;

% Particles receptors X
mesh = mshSphere(Nx,1);
X    = mesh.vtx;

% Particles transmitters Y (=X or not)
mesh = mshCube(Ny,2*[1 1 1]);
Y    = X;%mesh.vtx;
Ny   = Nx;%size(Y,1);

% Type
type = 'double';

% Wave number or frequency (Hz)
k = 5;
f = (k*340)/(2*pi);

% Green kernel -> exp(1i*k*r)/r
rxy   = @(X,Y) sqrt( (X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2 + (X(:,3)-Y(:,3)).^2 );
green = @(X,Y) exp(1i*k*rxy(X,Y))./(rxy(X,Y)+1e-8) .* (rxy(X,Y)>1e-8) + ...
    (sqrt(Nx)+1i*k) .* (rxy(X,Y)<=1e-8);

% Particles charges (multiples)
V    = (-1+2*rand(Ny,2)) + (-1+2i*rand(Ny,2));

% Spatial representation of particles
figure
plot3(X(:,1),X(:,2),X(:,3),'*b',Y(:,1),Y(:,2),Y(:,3),'*r')
axis equal 


%%% Full matrix computation
disp('~~~~~~~~~~~~~ EXACT MATRIX AS REFERENCE~~~~~~~~~~~~~')
tic
M = zeros(Nx,Ny,type);
for i = 1:Nx
    M(i,:) = green(X(i,:),Y).';
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


%%% H-Matrix computation
disp('~~~~~~~~~~~~~ H-MATRIX ~~~~~~~~~~~~~')
tic
Mh = hmx(X,Y,green,tol);
toc

tic
Mh2 = hmx(X,Y,M,tol);
toc

tic
Mh3 = hmx(Mh,0.1);
toc

tic
Ih = hmx(double(X),double(Y),I,tol);
toc

disp(' ')


%%% H-Matrix structure
disp('~~~~~~~~~~~~~ SPY STRUCTURE ~~~~~~~~~~~~~')
tic
figure
spy(Mh);
toc

tic
figure
spy(Mh2);
toc

tic
figure
spy(Mh3);
toc

tic
figure
spy(Ih);
toc

disp(' ')


%%% Full conversion
disp('~~~~~~~~~~~~~ FULL CONVERSION ~~~~~~~~~~~~~')
tic
sol = full(Mh);
toc
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = full(Mh2);
toc
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = full(Mh3);
toc
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = full(Ih);
toc
ref = full(I);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Sparse conversion
disp('~~~~~~~~~~~~~ SPARSE CONVERSION ~~~~~~~~~~~~~')
tic
sol = sparse(double(Mh));
toc
ref = double(M);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = sparse(Ih);
toc
ref = I;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Sparse conversion
disp('~~~~~~~~~~~~~ SPARSIFICATION ~~~~~~~~~~~~~')
tic
sol = sparse(double(Mh),sparse(M));
toc
ref = double(M);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = sparse(Ih,I);
toc
ref = I;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Low-rank conversion
disp('~~~~~~~~~~~~~ LOW-RANK CONVERSION ~~~~~~~~~~~~~')
tic
[A,B] = lowrank(Mh);
sol   = A*B;
toc
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Full conversion
disp('~~~~~~~~~~~~~ SUB-MATRIX ~~~~~~~~~~~~~')
idx = floor(Nx/4):floor(3*Nx/4);
jdx = floor(Ny/3):floor(2*Ny/3);

tic
sol = full(Mh,idx,jdx);
toc
ref = M(idx,jdx);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
[A,B] = lowrank(Mh,idx,jdx);
sol   = A*B; 
toc
ref = M(idx,jdx);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = sparse(Ih,idx,jdx);
toc
ref = I(idx,jdx);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
[A,B] = lowrank(Ih,idx,jdx);
sol   = A*B; 
toc
ref = I(idx,jdx);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Diagonal
disp('~~~~~~~~~~~~~ DIAGONAL ~~~~~~~~~~~~~')
tic
sol = diag(Mh);
toc
ref = diag(M);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = diag(Ih);
toc
ref = diag(I);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Single
disp('~~~~~~~~~~~~~ SINGLE ~~~~~~~~~~~~~')
tic
tmp = single(Mh);
toc
sol = full(tmp);
ref = single(M);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Double
disp('~~~~~~~~~~~~~ DOUBLE ~~~~~~~~~~~~~')
tic
tmp = double(Mh);
toc
sol = full(tmp);
ref = double(M);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Transposition
disp('~~~~~~~~~~~~~ TRANSPOSITION ~~~~~~~~~~~~~')
tic
tmp = Mh.';
toc
sol = full(tmp);
ref = M.';
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = Ih.';
toc
sol = full(tmp);
ref = I.';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Conjugate transposition
disp('~~~~~~~~~~~~~ CTRANSPOSITION ~~~~~~~~~~~~~')
tic
tmp = Mh';
toc
sol = full(tmp);
ref = M';
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = Ih';
toc
sol = full(tmp);
ref = I';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix right product
disp('~~~~~~~~~~~~~ H-MATRIX * MATRIX ~~~~~~~~~~~~~')
tic
sol = Mh * V;
toc
ref = MV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = Ih * double(V);
toc
ref = IV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = double(Mh) * I.';
toc
ref = double(M) * I.';
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = Ih * I.';
toc
ref = I * I.';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Left product H-Matrix
disp('~~~~~~~~~~~~~ MATRIX * H-MATRIX ~~~~~~~~~~~~~')
tmp = Mh.';
tic
sol = V.' * tmp;
toc
ref = MV.';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Scalar product
disp('~~~~~~~~~~~~~ SCALAR PRODUCT ~~~~~~~~~~~~~')
tic
tmp = (sqrt(2) .* Mh .* pi);
toc
sol = tmp * V;
ref = sqrt(2) * pi * MV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = (sqrt(2) .* Ih .* pi);
toc
sol = tmp * double(V);
ref = sqrt(2) * pi * IV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Uminus
disp('~~~~~~~~~~~~~ UMINUS ~~~~~~~~~~~~~')
tic
tmp = - Mh;
toc
sol = tmp * V;
ref = - MV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix addition
disp('~~~~~~~~~~~~~ H-MATRIX ADDITION/SUBSTRACTION ~~~~~~~~~~~~~')
tic
tmp = 2.*M + Mh - M;
toc
sol = tmp * V;
ref = 2 * MV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = 2.*I + Mh - I;
toc
sol = tmp * V;
ref = MV + IV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = 2.*M + Ih - M;
toc
sol = tmp * V;
ref = MV + IV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = 2.*I + Ih - I;
toc
sol = tmp * double(V);
ref = 2 * IV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix concatenation
disp('~~~~~~~~~~~~~ H-MATRIX CONCATENATION ~~~~~~~~~~~~~')
tic
tmp = [Mh Ih ; Ih Mh];
toc
sol = full(tmp);
tic
ref = [M I ; I M];
toc
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)


%%% H-Matrix product
disp('~~~~~~~~~~~~~ H-MATRIX PRODUCT ~~~~~~~~~~~~~')
tic
tmp = Mh * Mh.';
toc
sol = full(tmp);
ref = M * M.';
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)

tic
tmp = Ih.' * double(Mh);
toc
sol = full(tmp);
ref = I.' * double(M);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = double(Mh) * Ih.';
toc
sol = full(tmp);
ref = double(M) * I.';
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = Ih * Ih.';
toc
sol = sparse(tmp);
ref = I * I.';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix plus product
disp('~~~~~~~~~~~~~ H-MATRIX PLUS PRODUCT ~~~~~~~~~~~~~')
tic
tmp = plusmtimes(Ih,-pi,Mh,Mh.');
toc
sol = full(tmp);
ref = I - pi*(M * M.');
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix inversion
disp('~~~~~~~~~~~~~ INVERSION ~~~~~~~~~~~~~')
tic
tmp = inv(Mh);
toc
sol = tmp * V;
ref = M \ V;
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)

tic
tmp = inv(Ih);
toc
sol = tmp * V;
ref = I \ V;
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)

disp(' ')


%%% H-Matrix cholesky
disp('~~~~~~~~~~~~~ CHOLESKY FACTORISATION ~~~~~~~~~~~~~')
tic
Uh = chol(Ih);
toc
sol = sparse(Uh'*Uh);
ref = I;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix LU
disp('~~~~~~~~~~~~~ LU FACTORISATION ~~~~~~~~~~~~~')
tic
[Lh,Uh] = lu(Mh);
toc
sol = full(Lh*Uh);
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')
figure
subplot(1,2,1)
spy(Lh)
subplot(1,2,2)
spy(Uh)

tic
[Lh,Uh] = lu(Ih);
toc
sol = sparse(Lh*Uh);
ref = I;
norm(ref-sol,'inf')/norm(ref,'inf')
figure
subplot(1,2,1)
spy(Lh)
subplot(1,2,2)
spy(Uh)

disp(' ')


%%% H-Matrix \
disp('~~~~~~~~~~~~~ EXACT SOLVER ~~~~~~~~~~~~~')
tic
sol =  Mh \ V;
toc
tic
ref = M \ V;
toc
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol =  Ih \ double(V);
toc
tic
ref = I \ double(V);
toc
norm(ref-sol,'inf')/norm(ref,'inf')

[Lh,Uh] = lu(Mh);
tic
sol =  (Uh \ (Lh \ Mh)) * V;
toc
tic
ref = (M \ M ) * V;
toc
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')




disp('~~> Michto gypsilab !')



