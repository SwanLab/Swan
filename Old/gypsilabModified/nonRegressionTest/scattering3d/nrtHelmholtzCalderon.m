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
%|    #    |   FILE       : nrtCalderonHelmholtz.m                        |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Calderon Identity for integrals operators     |
%|  `---'  |                (poly Abboud Terrase p.194)                   |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Mesh sphere
mesh = mshSphere(1e3,1);

% Mesh cube
% mesh = mshCube(2e3,[1 1 1]);
% mesh = mesh.bnd

% Wave number
stp = mesh.stp;
k   = 1/stp(2);
f   = (k*340)/(2*pi);

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% grady Green kernel function --> G(x,y) = grady[exp(ik|x-y|)/|x-y|]
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Quadrature
sigma = dom(mesh,3);

% Finite element
u = fem(mesh,'P1');

% Mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,u);

% Single layer --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy
tic
S = 1/(4*pi) .* integral(sigma,sigma,u,Gxy,u);
S = S + 1/(4*pi) .* regularize(sigma,sigma,u,'[1/r]',u);
toc

% Double layer --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy
tic
D = 1/(4*pi) .* integral(sigma,sigma,u,dyGxy,ntimes(u));
D = D + 1/(4*pi) .* regularize(sigma,sigma,u,'grady[1/r]',ntimes(u));
toc

% Double layer transpose --> \int_Sx \int_Sy psi(x)' dnx G(x,y) psi(y) dx dy
Dt = D.';

% Hypersingular --> k^2 * \int_Sx \int_Sy n.psi(x) G(x,y) n.psi(y) dx dy
%                   - \int_Sx \int_Sy nxgrad(psi(x)) G(x,y) nxgrad(psi(y)) dx dy
tic
H = 1/(4*pi) .* (k^2 * integral(sigma,sigma,ntimes(u),Gxy,ntimes(u)) ...
    - integral(sigma,sigma,nxgrad(u),Gxy,nxgrad(u)));
H = H + 1/(4*pi) .* (k^2 * regularize(sigma,sigma,ntimes(u),'[1/r]',ntimes(u)) ...
    - regularize(sigma,sigma,nxgrad(u),'[1/r]',nxgrad(u)));
toc

% Calderon operators
Z = sparse(length(u),length(u));
I = [Id Z ; Z Id];
C = [-D S ; -H Dt];

% Calderon operator for interior problem (projector)
Ci = 0.5*I + C;
norm(Ci*(I\Ci)-Ci,'inf')/norm(Ci,'inf')

% Calderon operator for exterior problem (projector)
Ce = 0.5*I - C;
norm(Ce*(I\Ce)-Ce,'inf')/norm(Ce,'inf')

% Calderon Identity 1
A = D * (Id\S);
B = S * (Id\Dt);
norm(A-B,'inf')/norm(B,'inf')

% Calderon Identity 2
A = H * (Id\D);
B = Dt * (Id\H);
norm(A-B,'inf')/norm(B,'inf')

% Calderon Identity 3
A = D * (Id\D) - 0.25*Id;
B = S*(Id\H);
norm(A-B,'inf')/norm(B,'inf')

% Calderon Identity 4
A = Dt * (Id\Dt) - 0.25*Id;
B = H*(Id\S);
norm(A-B,'inf')/norm(B,'inf')





disp('~~> Michto gypsilab !')



