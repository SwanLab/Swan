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
%|    #    |   FILE       : nrtHmxCompressorTotal.m                       |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Compare compressor with total pivoting        |
%|  `---'  |                                                              |
%+========================================================================+

clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Precision
tol = 1e-2

% Read
M   = imread('chevre.jpg');
M   = double(M(:,:,1));
Msp = sparse(M);

% Graphical representation
figure
imagesc(M)
axis equal off
colormap gray

figure
spy(Msp)


%%% FULL COMPRESSION
% Compression SVD
disp('------> FULL SVD')
tic
[A,B] = hmxSVD(M,tol);
toc
size(A)
norm(M-A*B,'fro')/norm(M,'fro')

% Compression RSVD
disp('------> FULL RSVD')
tic
rk = rank(M)
toc
tic
[A,B] = hmxRSVD(M,tol,rk);
toc
size(A)
norm(M-A*B,'fro')/norm(M,'fro')

% Compression ACA
disp('------> FULL ACA')
tic
[A,B] = hmxACA(M,tol);
toc
size(A)
norm(M-A*B,'fro')/norm(M,'fro')

% Compression RSVD low-rank
disp('------> FULL RSVD')
tic
[Ap,Bp] = hmxRSVD(A,B,tol);
toc
size(Ap)
norm(M-Ap*Bp,'fro')/norm(M,'fro')

% Compression QRSVD
disp('------> FULL QRSVD')
tic
[A,B] = hmxQRSVD(A,B,tol);
toc
size(A)
norm(M-A*B,'fro')/norm(M,'fro')

% Graphical representation
figure
imagesc(A*B)
axis equal off
colormap gray


%%% SPARSE COMPRESSION
% Compression RSVD
disp('------> SPARSE RSVD')
tic
rk = sprank(Msp)
toc
[A,B] = hmxRSVD(Msp,tol,rk);
toc
size(A)
norm(M-A*B,'fro')/norm(M,'fro')

% Compression ACA
disp('------> SPARSE ACA')
tic
[A,B] = hmxACA(Msp,tol);
toc
size(A)
% norm(M-A*B,'fro')/norm(M,'fro')   % to be debbugged

% Compression RSVD low-rank
disp('------> SPARSE RSVD')
tic
[Ap,Bp] = hmxRSVD(A,B,tol);
toc
size(Ap)
norm(M-Ap*Bp,'fro')/norm(M,'fro')

% Compression QRSVD
disp('------> SPARSE QRSVD')
tic
[Asp2,Bsp2] = hmxQRSVD(A,B,tol);
toc
size(Asp2)
norm(M-Asp2*Bsp2,'fro')/norm(M,'fro')

% Graphical representation
figure
imagesc(A*B)
axis equal off
colormap gray



disp('~~> Michto gypsilab !')
