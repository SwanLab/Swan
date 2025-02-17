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
%|    #    |   FILE       : nrtFfmHelmholtzBWdir.m                        |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Solve dirichlet scatering problem with        |
%|  `---'  |                Brackage-Werner formulation                   |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N   = 1e3
tol = 1e-3
typ = 'P1'
gss = 3
X0  = [0 0 -1]

% Spherical mesh
sphere = mshSphere(N,1);
sigma  = dom(sphere,gss);    
figure
plot(sphere)
axis equal

% Radiative mesh
square     = mshSquare(5*N,[5 5]);
square.vtx = [square.vtx(:,1) zeros(size(square.vtx,1),1) square.vtx(:,2)];
hold on
plot(square)

% Frequency adjusted to maximum esge size
stp = sphere.stp;
k   = 1/stp(2)
f   = (k*340)/(2*pi)

% Incident wave
PW = @(X) exp(1i*k*X*X0');

% Incident wave representation
plot(sphere,real(PW(sphere.vtx)))
plot(square,real(PW(square.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
view(0,10)


%%% PREPARE OPERATORS
disp('~~~~~~~~~~~~~ PREPARE OPERATORS ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy      = '[exp(ikr)/r]';
gradyGxy = {'grady[exp(ikr)/r]1','grady[exp(ikr)/r]2','grady[exp(ikr)/r]3'};

% Finite elements
u = fem(sphere,typ);
v = fem(sphere,typ);

% Coupling coeff
beta = 1i*k*0.5;

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,v);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sr = 1/(4*pi) .* regularize(sigma,sigma,u,'[1/r]',v);
S  = 1/(4*pi) .* integral(sigma,sigma,u,Gxy,k,v,tol) + Sr;
toc

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy 
tic
Dr = 1/(4*pi) .* regularize(sigma,sigma,u,'grady[1/r]',ntimes(v));
D  = 1/(4*pi) .* integral(sigma,sigma,u,gradyGxy,k,ntimes(v),tol) + Dr;
toc

% Final operator : [1i*k*beta*S - (Id/2 + D)]
tic
LHS  = beta.*S - (0.5*Id + D);
toc

% Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
RHS = - integral(sigma,u,PW);


%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% Preconditionneur ILU
tic
[L,U] = ilu(beta.*Sr - (0.5*Id + Dr));
toc

% Solve linear system : [1i*k*beta*S - (Id/2 + D)] = P0
tic
mu = mgcr(@(V) LHS*V,RHS,[],tol,100,L,U);
toc

% Jump for derivative
lambda = beta * mu;


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
theta = 2*pi/1e3 .* (1:1e3)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
Ginf      = '[exp(-ikxy)]';
gradxGinf = {'gradx[exp(-ikxy)]1','gradx[exp(-ikxy)]2','gradx[exp(-ikxy)]3'};

% Finite element infinite operators
Sinf = 1/(4*pi) .* integral(nu,sigma,Ginf,k,v,tol);
Dinf = 1/(4*pi) .* integral(nu,sigma,gradxGinf,k,ntimes(v),tol);

% Finite element radiation  
sol = Sinf*lambda - Dinf*mu;

% Analytical solution
ref = sphereHelmholtz('inf','dir',1,k,nu); 
norm(ref-sol,2)/norm(ref,2)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation
figure
plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')


%%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,v);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sbnd = 1/(4*pi) .* (integral(sigma,sigma,u,Gxy,k,v,tol) + ...
    regularize(sigma,sigma,u,'[1/r]',v));
toc

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy 
tic
Dbnd = 1/(4*pi) .* (integral(sigma,sigma,u,gradyGxy,k,ntimes(v),tol) + ...
    regularize(sigma,sigma,u,'grady[1/r]',ntimes(v)));
toc

% Boundary solution
Psca = Id\(Sbnd*lambda - (0.5*Id*mu + Dbnd*mu));
Pinc = PW(u.dof);
Pbnd = Pinc + Psca;

% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy 
tic
Sdom = 1/(4*pi) .* (integral(square.vtx,sigma,Gxy,k,v,tol) + ...
    regularize(square.vtx,sigma,'[1/r]',v));
toc

% Finite element radiative operator --> \int_Sx \int_Sy psi(x)' grady(G(x,y)) ny.psi(y) dx dy 
tic
Ddom = 1/(4*pi) .* ( integral(square.vtx,sigma,gradyGxy,k,ntimes(v),tol) + ...
    regularize(square.vtx,sigma,'grady[1/r]',ntimes(v)) );
toc

% Domain solution
Psca = Sdom*lambda - Ddom*mu;
Pinc = PW(square.vtx);
Pdom = Pinc + Psca;

% Annulation sphere interieure
r             = sqrt(sum(square.vtx.^2,2));
Pdom(r<=1.01) = Pinc(r<=1.01);

% Graphical representation
figure
plot(sphere,abs(Pbnd))
axis equal;
hold on
plot(square,abs(Pdom))
title('Total field solution')
colorbar
hold off
view(0,10)


%%% ANAYTICAL SOLUTIONS FOR COMPARISONS
% Analytical solution
Pbnd = sphereHelmholtz('dom','dir',1,k,1.001*sphere.vtx) + PW(sphere.vtx);
Pdom = sphereHelmholtz('dom','dir',1,k,square.vtx) + PW(square.vtx);

% Solution representation
figure
plot(sphere,abs(Pbnd))
axis equal;
hold on
plot(square,abs(Pdom))
title('Analytical solution')
colorbar
hold off
view(0,10)



disp('~~> Michto gypsilab !')


