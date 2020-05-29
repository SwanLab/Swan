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
%|    #    |   FILE       : nrtHmxLowrank.m                               |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : Compare compressor with total pivoting        |
%|  `---'  |                                                              |
%+========================================================================+

clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Dimensions
Nx  = 100;
Ny  = 101;
rk  = 10;
tol = 1e-3;

% Compress zeros
M     = zeros(Nx,Ny);
[A,B] = hmxACA(M,tol);
norm(A*B,'inf')

% Compress ones
M     = ones(Nx,Ny);
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Compress eye 
M          = eye(Nx,Ny);
[~,~,flag] = hmxACA(M,tol);
flag

% Compress circular permutation for eye
M = circshift(eye(Nx,Ny),3,2);
[~,~,flag] = hmxACA(M,tol);
flag

% Compress 1-eye 
M          = 1-eye(Nx,Ny);
[~,~,flag] = hmxACA(M,tol);
flag

% Compress random
M     = rand(Nx,rk) * rand(rk,Ny);
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Compress random complex
M     = (-1+2*rand(Nx,rk)+1i*(-1+2*rand(Nx,rk))) * (rand(rk,Ny) + 1i*rand(rk,Ny));
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Compress 1x2 block with zeros
M     = rand(Nx,rk) * rand(rk,Ny);
Z     = zeros(Nx,10);
M     = [M Z];
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Compress 2x1 block with zeros
M     = rand(Nx,rk) * rand(rk,Ny);
Z     = zeros(10,Ny);
M     = [M ; Z];
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Compress 2x2 block random
M     = rand(Nx,rk) * rand(rk,Ny);
Z     = zeros(Nx,Ny);
M     = [M Z ; Z M];
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Compress 4x4 block random
M     = rand(Nx,rk) * rand(rk,Ny);
Z     = zeros(Nx,Ny);
M     = [M Z ; Z M];
Z     = zeros(size(M));
M     = [M Z ; Z M];
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Compress 2x2 inegal block random (to be hardly improved...)
M     = rand(Nx,rk) * rand(rk,Ny);
Z     = zeros(Nx,1);
M     = [M zeros(Nx,1) ; zeros(1,Ny) 1];
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Compress 4x4 inegal block random (to be hardly improved...)
M     = rand(Nx,rk) * rand(rk,Ny);
Z     = zeros(Nx,1);
M     = [M zeros(Nx,3) ; zeros(3,Ny) eye(3)];
[A,B] = hmxACA(M,tol);
norm(A*B-M,'inf')/norm(M,'inf')

% Build orthogonal mesh
vtx   = [0 0 0 ; eye(3)];
elt   = [1 3 2 ; 1 2 4];
col   = [1 ; 2];
mesh  = msh(vtx,elt,col);
mesh  = refine(mesh,0.2);
[~,I] = sort(mesh.col);
mesh  = mesh.sub(I);

% Translate
mesh2          = mesh;
mesh2.vtx(:,1) = mesh2.vtx(:,1) + 2;

% Graphical representation
figure
plot(mesh)
hold on
plotNrm(mesh)
plot(mesh2)
axis equal
view(-180,30)

% Build double layer matrix
Gxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',1);
Gxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',1);
Gxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',1);
D      = integral(dom(mesh,3),dom(mesh2,3),fem(mesh,'P0'),Gxy,ntimes(fem(mesh2,'P0')));

% Graphical represenation
figure
imagesc(abs(D))
colorbar

% Compression
[A,B] = hmxACA(D,tol);
norm(A*B-D,'inf')/norm(D,'inf')




disp('~~> Michto gypsilab !')




