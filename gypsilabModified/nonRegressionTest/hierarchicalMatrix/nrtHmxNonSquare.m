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
%|    #    |   FILE       : nrtHmxNonSquare.m                             |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Marc Bakry                  |
%|  ( # )  |   CREATION   : 14.01.2019                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : H-Matrix algebra with non square matrix       |
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
Ny = 1e2;

% Accuracy
tol = 1e-6;

% Particles receptors X
mesh = mshSphere(Nx,1);
X    = mesh.vtx;

% Particles transmitters Y (=X or not)
mesh = mshSphere(Ny,1);
Y    = mesh.vtx;

% Wave number or frequency (Hz)
k = 5;
f = (k*340)/(2*pi);

% Green kernel -> exp(1i*k*r)/r
green = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Particles charges (multiples)
Vx = (-1+2*rand(Nx,2)) + (-1+2i*rand(Nx,2));
Vy = (-1+2*rand(Ny,2)) + (-1+2i*rand(Ny,2));

% Spatial representation of particles
figure
plot3(X(:,1),X(:,2),X(:,3),'*b',Y(:,1),Y(:,2),Y(:,3),'*r')
axis equal 

% Full matrix
tic
Mxx = zeros(Nx,Nx);
Mxy = zeros(Nx,Ny);
Myy = zeros(Ny,Ny);
for i = 1:Nx
    Mxx(i,:) = green(X(i,:),X).';
    Mxy(i,:) = green(X(i,:),Y).';
end
for i = 1:Ny
    Myy(i,:) = green(Y(i,:),Y).';
end
toc

% H-Matrix
tic
Mxxh = hmx(X,X,green,tol);
toc
figure
spy(Mxxh);

tic
Mxyh = hmx(X,Y,green,tol);
toc
figure
spy(Mxyh);

tic
Myyh = hmx(Y,Y,green,tol);
toc
figure
spy(Myyh);


% H-Matrix-vector product
ref = Mxx * Vx;
sol = Mxxh * Vx;
norm(ref-sol,'inf')/norm(ref,'inf')

ref = Mxy * Vy;
sol = Mxyh * Vy;
norm(ref-sol,'inf')/norm(ref,'inf')

ref = Myy * Vy;
sol = Myyh * Vy;
norm(ref-sol,'inf')/norm(ref,'inf')


% H-Matrix product
ref = (Mxx * Mxy * Myy) * Vy;
sol = (Mxxh * Mxyh * Myyh) * Vy;
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(Mxxh * Mxyh * Myyh)

ref = (Myy * Mxy.' * Mxx) * Vx;
sol = (Myyh * Mxyh.' * Mxxh) * Vx;
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(Myyh * Mxyh.' * Mxxh)

ref = (Mxy * Mxy.') * Vx;
sol = (Mxyh * Mxyh.') * Vx;
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(Mxyh * Mxyh.')


% H-Matrix LU
ref = Mxx \ Mxy;
sol = Mxxh \ Mxyh;
norm(ref*Vy-sol*Vy,'inf')/norm(ref*Vy,'inf')

ref = (Mxy.') / (Mxx);
sol = (Mxyh.') / (Mxxh);
norm(ref*Vx-sol*Vx,'inf')/norm(ref*Vx,'inf')



disp('~~> Michto gypsilab !')



