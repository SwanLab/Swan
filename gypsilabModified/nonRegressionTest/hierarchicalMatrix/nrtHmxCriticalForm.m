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
%|    #    |   FILE       : nrtHmxCriticalForm.m                          |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 01.05.2019                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : Build particles H-Matrix for critical geometry|
%|  `---'  |                                                              |
%+========================================================================+

% Nettoyage
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Dimensions
N = 1000;

% Accuracy
tol = 1e-3;

% Meshes
mesh = mshSquare(N,[1 1]);
N    = size(mesh.vtx,1);

% Particles
X      = mesh.vtx;
Y      = mesh.vtx;
Y(:,1) = Y(:,1) + 1;
Y(:,2) = Y(:,2) + 0;
Y(:,3) = Y(:,3) + 0.2;

% Wave number or frequency (Hz)
k = 5;
f = (k*340)/(2*pi);

% Green kernel -> exp(1i*k*r)/r
rxy   = @(X,Y) sqrt( (X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2 + (X(:,3)-Y(:,3)).^2 );
green = @(X,Y) exp(1i*k*rxy(X,Y))./(rxy(X,Y)+1e-8) .* (rxy(X,Y)>1e-8) + ...
    (sqrt(N)+1i*k) .* (rxy(X,Y)<=1e-8);

% Charges
V = -1+2*rand(N,1);

% Spatial representation of particles
figure
plot(msh(X),'b')
hold on
plot(msh(Y),'r')
view(45,45)
axis equal 

% Full product
tic
ref = zeros(N,1);
for i = 1:N
    ref(i) = green(X(i,:),Y).' * V;
end
toc

% Admissibility
isfar = hmxFar(hmx(X,Y,tol))

% Compression
[A,B] = hmxACA(X,Y,green,tol);
size(A)
sol = A * (B * V);
norm(ref-sol)/norm(ref)

% H-Matrix
Mh  = hmx(X,Y,green,tol);
sol = Mh * V;
norm(ref-sol)/norm(ref)

% Graphical representation
figure
spy(Mh)





disp('~~> Michto gypsilab !')





