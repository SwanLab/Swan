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
%|    #    |   FILE       : nrtHmxCompressorPartial.m                     |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Compare compressor with partial pivoting      |
%|  `---'  |                                                              |
%+========================================================================+

clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Fibonacci sphere 
N    = 1e3;
mesh = mshSphere(N,1);
X    = mesh.vtx;

% Potentiel aleatoire aux emmeteurs
V  = (-1-1i) + (2+2i)*rand(N,1);

% Representation graphique
figure(1)
plot3(X(:,1),X(:,2),X(:,3),'*b')
grid on
axis equal

% Noyau de Fourier
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3);
green = @(X,Y,k) exp(1i*k*xdoty(X,Y));

% Initialisation
K   = 1:5;
rnk = zeros(6,length(K));

% Boucle sur k
for i = 1:length(K)
    % Nombre d'onde
    k = K(i);
   
    % Full green kernel
    tic
    [I,J] = ndgrid(1:N,1:N);
    Gxy   = green(X(I,:),X(J,:),k);
    Gxy   = reshape(Gxy,N,N);
    toc

    % SVD
    [A,B] = hmxSVD(Gxy,1e-3);
    toc
    size(A)
    rnk(1,i) = size(A,2);
    norm(A*B-Gxy,'fro')./norm(Gxy,'fro')
    
    % Compression SVD
    tic
    [A,B] = hmxSVD(Gxy,1e-6);
    toc
    size(A)
    rnk(2,i) = size(A,2);
    norm(A*B-Gxy,'fro')./norm(Gxy,'fro')
    
    % Compression SVD
    tic
    [A,B] = hmxSVD(Gxy,1e-9);
    toc
    rnk(3,i) = size(A,2);
    norm(A*B-Gxy,'fro')./norm(Gxy,'fro')    

    % Compression ACA, pivotage partiel
    tic
    [A,B] = hmxACA(X,X,@(X,Y) green(X,Y,k),1e-3);
    toc
    size(A)
    rnk(4,i) = size(A,2);
    norm(A*B-Gxy,'fro')./norm(Gxy,'fro')
    
    % Recompression RSVD
    tic
    [Ap,Bp] = hmxRSVD(A,B,1e-3);
    toc
    size(Ap)
    norm(Ap*Bp-Gxy,'fro')./norm(Gxy,'fro')
    
    % Recompression QRSVD
    tic
    [Ap,Bp] = hmxQRSVD(A,B,1e-3);
    toc
    size(Ap)
    norm(Ap*Bp-Gxy,'fro')./norm(Gxy,'fro')  
    
    % Compression ACA, pivotage partiel
    tic
    [A,B] = hmxACA(X,X,@(X,Y) green(X,Y,k),1e-6);
    toc
    rnk(5,i) = size(A,2);
    norm(A*B-Gxy,'fro')./norm(Gxy,'fro')
    
    % Compression ACA, pivotage partiel
    tic
    [A,B] = hmxACA(X,X,@(X,Y) green(X,Y,k),1e-9);
    toc
    rnk(6,i) = size(A,2);
    norm(A*B-Gxy,'fro')./norm(Gxy,'fro')
    '====================================='
end

% Rang en fonction de k
figure(3)
plot(K,rnk,'+-')
legend({'svd at 1e-3','svd at 1e-6','svd at 1e-9',...
    'aca at 1e-3','aca at 1e-6','aca at 1e-9'})
xlabel('Wave number k')
ylabel('Rank')
grid on



disp('~~> Michto gypsilab !')

