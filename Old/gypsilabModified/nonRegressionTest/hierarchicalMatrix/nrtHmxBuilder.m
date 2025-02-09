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
%|    #    |   FILE       : nrtHmxBuilder.m                               |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Build H-Matrix and compare to full product    |
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

% Wave number or frequency (Hz)
k = 5;
f = (k*340)/(2*pi);

% Green kernel -> exp(1i*k*r)/r
rxy   = @(X,Y) sqrt( (X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2 + (X(:,3)-Y(:,3)).^2 );
green = @(X,Y) exp(1i*k*rxy(X,Y))./(rxy(X,Y)+1e-8) .* (rxy(X,Y)>1e-8) + ...
    (sqrt(Nx)+1i*k) .* (rxy(X,Y)<=1e-8);
% green = @(X,Y) femGreenKernel(X,Y,'[1/r]',k);

% Particles charges
V = (-1+2*rand(Ny,1)) + (-1+2i*rand(Ny,1));

% Spatial representation of particles
figure
plot3(X(:,1),X(:,2),X(:,3),'*b',Y(:,1),Y(:,2),Y(:,3),'*r')
axis equal 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCUL DIRECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FULL PRODUCT ~~~~~~~~~~~~~')
tic

% Recepteurs au hasard
Nt   = 100;
ind  = 1:ceil(Nx/Nt):Nx;
Xloc = X(ind,:);

% Produit Matrice-Vecteur sur tous les emeteurs
ref = zeros(length(ind),1);
for i = 1:length(ind)
    ref(i) = green(Xloc(i,:),Y).' * V;
end
toc
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FAST PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FAST PRODUCT ~~~~~~~~~~~~~')

% Hierarchical Matrix
tic
Mh = hmx(X,Y,green,tol);
toc

% Graphical representation
figure
subplot(2,2,1:2)
spy(Mh)

% Matrix vector product
tic
sol = Mh * V;
toc

% Erreur
norm(ref-sol(ind))/norm(ref)
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LU FACTORIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ LU FACTORIZATION ~~~~~~~~~~~~~')

% Factorization LU
tic
[Lh,Uh] = lu(Mh);
toc

% Graphical representation
subplot(2,2,3); spy(Lh)
subplot(2,2,4); spy(Uh)

% Matrix vector product
sol = Lh * (Uh * V);

% Error
norm(ref-sol(ind))/norm(ref)
disp(' ')


disp('~~> Michto gypsilab !')
