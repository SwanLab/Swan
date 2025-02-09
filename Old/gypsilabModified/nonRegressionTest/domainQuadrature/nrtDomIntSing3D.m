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
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtDomIntSing3D.m                             |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 05.09.2018                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Singular integration over a triangle          |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Reference mesh (triangle)
vtx  = 1.0 * [0 0 0 ; 1 0 0 ; 0 1 0];
elt  = [1 2 3];
mesh = msh(vtx,elt);
edge = mesh.edg;

% Observation points (100 random & particular points)
N = 100;
X = [-2 + 5*rand(N,3);
    mesh.vtx;
    mesh.ctr
    -mesh.ctr
    edge.ctr;
    -edge.ctr];

% Graphical representation
plot(mesh,1)
hold on
plot3(X(:,1),X(:,2),X(:,3),'or')
axis equal
view(30,30)

% Analytical integration (exact)
tic
S   = mesh.vtx;
n   = mesh.nrm;
tau = cell2mat(mesh.tgt');
nu  = cell2mat(mesh.nrmEdg');
[Rm1a,rRm1a,gradRm1a,gradrRm1a] = domSemiAnalyticInt3D(X,S,n,tau,nu);
toc
   
% 2D simpson integration ("exact")
tic
Rm1s      = zeros(size(Rm1a));
rRm1s     = zeros(size(rRm1a));
gradRm1s  = zeros(size(gradRm1a));
gradrRm1s = zeros(size(gradrRm1a));
for i = 1:size(X,1)
    % Function |r| in the 2D plane of the triangle
    fun = @(x1,x2) sqrt( (x1-X(i,1)).^2 + (x2-X(i,2)).^2 + X(i,3).^2 );
    
    % Scalar integration 1/|r|
    Rm1s(i) = integral2(@(x1,x2) 1./fun(x1,x2), ...
        0, 1 , 0, @(x1) 1-x1);
    
    % Vectorial integration r/|r|
    rRm1s(i,1) = integral2(@(x1,x2) (x1-X(i,1))./fun(x1,x2), ...
        0, 1 , 0, @(x1) 1-x1);
    rRm1s(i,2) = integral2(@(x1,x2) (x2-X(i,2))./fun(x1,x2), ...
        0, 1 , 0, @(x1) 1-x1);
    rRm1s(i,3) = -X(i,3) .*  Rm1s(i);
    
    % Vectorial integration grad(1/|r|) = -r/|r|^3
    if (X(i,3) == 0)
        gradRm1s(i,:) = 0;
    else
        gradRm1s(i,1) = integral2(@(x1,x2) - (x1-X(i,1))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradRm1s(i,2) = integral2(@(x1,x2) - (x2-X(i,2))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradRm1s(i,3) = integral2(@(x1,x2) X(i,3)./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
    end
    
    % Tensorial integration grad(r/|r|) = Id/r - rxr/|r|^3
    if (X(i,3) == 0)
        gradrRm1s(i,:) = 0;
    else
        gradrRm1s(i,1,1) = integral2(@(x1,x2) 1./fun(x1,x2) - (x1-X(i,1)).*(x1-X(i,1))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradrRm1s(i,2,1) = integral2(@(x1,x2) - (x2-X(i,2)).*(x1-X(i,1))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradrRm1s(i,3,1) = integral2(@(x1,x2) - (0-X(i,3)).*(x1-X(i,1))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradrRm1s(i,1,2) = integral2(@(x1,x2) - (x1-X(i,1)).*(x2-X(i,2))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradrRm1s(i,2,2) = integral2(@(x1,x2) 1./fun(x1,x2) - (x2-X(i,2)).*(x2-X(i,2))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradrRm1s(i,3,2) = integral2(@(x1,x2) - (0-X(i,3)).*(x2-X(i,2))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradrRm1s(i,1,3) = integral2(@(x1,x2) - (x1-X(i,1)).*(0-X(i,3))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradrRm1s(i,2,3) = integral2(@(x1,x2) - (x2-X(i,2)).*(0-X(i,3))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
        gradrRm1s(i,3,3) = integral2(@(x1,x2) 1./fun(x1,x2) - (0-X(i,3)).*(0-X(i,3))./fun(x1,x2).^3, ...
            0, 1 , 0, @(x1) 1-x1);
    end
end
toc

% Relative errors 1/|r|
figure
semilogy( abs(Rm1a-Rm1s)./abs(Rm1a) )
title('Relative error (log) : \int 1/|r|')
grid on
I = 1:N;
disp('Relative error (Linf) : \int 1/|r|')
norm(Rm1a(I)-Rm1s(I),'inf')./norm(Rm1a(I),'inf')
norm(Rm1a-Rm1s,'inf')./norm(Rm1a,'inf')

% Relative errors r/|r|
figure
semilogy( abs(rRm1a-rRm1s)./abs(rRm1a) )
title('Relative error (log) : \int r/|r|')
grid on
disp('Relative error (Linf) : \int r/|r|')
norm(rRm1a(I,:)-rRm1s(I,:),'inf')./norm(rRm1a(I,:),'inf')

% Relative errors grad(1/|r|)
figure
semilogy( abs(gradRm1a-gradRm1s)./abs(gradRm1a) )
title('Relative error (log) : \int grad(1/|r|)')
grid on
disp('Relative error (Linf) : \int grad(1/|r|)')
norm(gradRm1a(I,:)-gradRm1s(I,:),'inf')./norm(gradRm1a(I,:),'inf')

% Relative errors gradx(r/|r|)
figure
semilogy( abs(gradrRm1a(:,:,1)-gradrRm1s(:,:,1))./abs(gradrRm1a(:,:,1)) )
title('Relative error (log) : \int grad(rx/|r|)')
grid on
disp('Relative error (Linf) : \int grad(rx/|r|)')
norm(gradrRm1a(I,:,1)-gradrRm1s(I,:,1),'inf')./norm(gradrRm1a(I,:,1),'inf')
disp('Relative error (Linf) : \int grad(ry/|r|)')
norm(gradrRm1a(I,:,2)-gradrRm1s(I,:,2),'inf')./norm(gradrRm1a(I,:,2),'inf')
disp('Relative error (Linf) : \int grad(rz/|r|)')
norm(gradrRm1a(I,:,3)-gradrRm1s(I,:,3),'inf')./norm(gradrRm1a(I,:,3),'inf')




disp('~~> Michto gypsilab !')




