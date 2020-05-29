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
%|    #    |   FILE       : nrtFfmBuilder.m                               |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Non regression test for convolution using     |
%|  `---'  |                arbitrary kenel with X and Y random sample    |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Unitary sphere
unit = @(u) u ./ sqrt(sum(u.^2,2)*ones(1,size(u,2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ DEFINITIONS ~~~~~~~~~~~~~')

% Data type
type = 'double'

% Green kernel
green = '[exp(ikr)/r]'

% Accuracy
tol = 1e-3

% Receivers positions
Nx = 1e3
X  = -1 + 2*rand(Nx,3,type);
X  = unit(X);

% Emmiters positions
Ny = 2e3
Y  = -1 + 2*rand(Ny,3,type);

% Wave number
k = 5

% Potential 
V = -(1+1i) + 2*(rand(Ny,2,type) + 1i*rand(Ny,2,type));

% Graphical representation
figure
plot3(X(:,1),X(:,2),X(:,3),'*b',Y(:,1),Y(:,2),Y(:,3),'*r')
axis equal 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FULL PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FULL PRODUCT ~~~~~~~~~~~~~')
tic

% Random receivers indices
Nt   = 100;
ind  = 1:ceil(Nx/Nt):Nx;
Xloc = X(ind,:);

% Full and exact matrix-vector product on randomized indices
ref = zeros(length(ind),size(V,2),type);
for i = 1:length(ind)
    ref(i,:) = ffmGreenKernel(Xloc(i,:),Y,green,k).' * V;
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FFM PRODUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FFM PRODUCT ~~~~~~~~~~~~~')

% Fast and free memory
tic
Mv = ffm(X,Y,green,k,tol);
toc

% Matrix vector product
tic
sol = Mv * V;
toc

% Compare to randomized indices values
norm(ref-sol(ind,:))/norm(ref)



disp('~~> Michto gypsilab !')



