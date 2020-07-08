%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as publiduP1ed by the Free Software Foundation (version 3 or  |
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
%|    #    |   FILE       : nrtFemOperators.m                             |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Francois Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Finites elements operators validation         |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parametres
N = 1e3;
L = [1 2 3];
V = prod(L);
S = 2*(L(1)*L(2) + L(1)*L(3) + L(2)*L(3));

% Meshes
mesh  = mshCube(N,L);
bound = mesh.bnd;

% Graphical representation
figure
plot(mesh)
alpha(0.1)
axis equal

figure
plot(bound)
hold on
plotNrm( bound)
alpha(0.5)
axis equal

% Domains
omega = dom(mesh,4);
sigma = dom(bound,3);

% Volumic finite elements
uP1   = fem(mesh,'P1');
uNED  = fem(mesh,'NED');
duP1  = fem(bound,'P1');
duNED = fem(bound,'NED');
duRWG = fem(bound,'RWG');


% \int_omega uP1 uP1
I = integral(omega,uP1,uP1);
error1 = sum(sum(I)) - V

% \int_sigma trace(uP1) trace(uP1)
I = integral(sigma,uP1,uP1);
error2 = sum(sum(I)) - S

% \int_sigma duP1 duP1 
I = integral(sigma,duP1,duP1);
error3 = sum(sum(I)) - S

% \int_omega grad(uP1) grad(uP1)
I  = integral(omega,grad(uP1),grad(uP1));
sum(sum(I))
I2 = integral(omega,grad(uP1,1),grad(uP1,1)) + ...
    integral(omega,grad(uP1,2),grad(uP1,2)) + ...
    integral(omega,grad(uP1,3),grad(uP1,3));
error4 = norm(I-I2,'inf')

% \int_sigma trace(grad(uP1)) trace(grad(uP1))
I = integral(sigma,grad(uP1),grad(uP1));
error5 = sum(sum(I))

% \int_sigma grad(duP1) grad(duP1)
I = integral(sigma,grad(duP1),grad(duP1));
error6 = sum(sum(I))

% \int_sigma grad(duP1) n*(duP1)
I = integral(sigma,grad(duP1),ntimes(duP1))
I = integral(sigma,grad(duP1,1),ntimes(duP1,1)) + ...
    integral(sigma,grad(duP1,2),ntimes(duP1,2)) + ...
    integral(sigma,grad(duP1,3),ntimes(duP1,3)) 

% \int_sigma duP1 n*(duP1)
N = sigma.qudNrm;
F = {@(X) N(:,1) , @(X) N(:,2) , @(X) N(:,3)};
A = integral(sigma,duP1,ntimes(duP1));
B = integral(sigma,duP1,F,duP1);
errorNrm1 = norm(A{1}-B{1},'inf')
errorNrm2 = norm(A{2}-B{2},'inf')
errorNrm3 = norm(A{3}-B{3},'inf')
C = integral(sigma,duP1,ntimes(duP1,1));
errorNrm1 = norm(A{1}-C,'inf')
C = integral(sigma,duP1,ntimes(duP1,2));
errorNrm2 = norm(A{2}-C,'inf')
C = integral(sigma,duP1,ntimes(duP1,3));
errorNrm3 = norm(A{3}-C,'inf')

% \int sigma duRWG nx(duRWG)
I = integral(sigma,duRWG,nx(duRWG));
sum(sum(I))

% \int n*(duP1) duRWG
I = integral(sigma,ntimes(duP1),duRWG)

% \int_sigma duP1 div(duRWG) = - \int_sigma grad(duP1) duRWG
A = integral(sigma,duP1,div(duRWG)); 
B = integral(sigma,grad(duP1),duRWG);
errorP1_divRWG = norm(A+B,'inf')

% \int_sigma duNED nxgrad(duP1) = \int_sigma rot(duNED) duP1
A = integral(sigma, duNED, nxgrad(duP1));
B = integral(sigma, curl(duNED), duP1);
errorNED_curlP1 = norm(A-B,'inf')

% \int_omega uNED rot(uNED) = \int_omega rot(u_NED) uNED - \int_sigma 0
uNED = dirichlet(uNED,bound);
A    = integral(omega,uNED,curl(uNED));
B    = integral(omega,curl(uNED),uNED);
C    = integral(sigma,nx(uNED),uNED)
norm(A-B,'inf')
uNED = dirichlet(uNED,[]);

% \int_sigma nx(duNED) nx(duNED) = \int_sigma duRWG duRWG 
A = integral(sigma,nx(duNED),duNED);
B = integral(sigma,duRWG,duNED);
norm(A-B,'inf')

% \int_sigma nx(uNED) uNED = P \int_sigma duRWG duNED P' 
A = integral(sigma,nx(uNED),uNED);
P = restriction(uNED,bound);
B = P' * integral(sigma,nx(duNED),duNED) * P;
norm(A-B,'inf')

% \int_omega uNED rot(uNED) = \int_omega rot(u_NED) uNED - \int_sigma nx(uNED) uNED
A = integral(omega,uNED,curl(uNED));
B = integral(omega,curl(uNED),uNED);
C = P' * integral(sigma,nx(duNED),duNED) * P;
D = integral(sigma,nx(uNED),uNED);
errorintparts = norm(A-B+C,'inf')
errorintparts = norm(A-B+D,'inf')



disp('~~> Michto gypsilab !')


