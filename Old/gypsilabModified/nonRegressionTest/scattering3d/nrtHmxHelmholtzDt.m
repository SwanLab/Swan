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
%|    #    |   FILE       : nrtHmxHelmholtzDt.m                           |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Solve neumann scatering problem with double   |
%|  `---'  |                layer transpose potential                     |
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
PW         = @(X) exp(1i*k*X*X0');
gradxPW{1} = @(X) 1i*k*X0(1) .* PW(X);
gradxPW{2} = @(X) 1i*k*X0(2) .* PW(X);
gradxPW{3} = @(X) 1i*k*X0(3) .* PW(X);

% Incident wave representation
plot(sphere,real(PW(sphere.vtx)))
plot(square,real(PW(square.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
view(0,10)


%%% PREPARE OPERATOR
disp('~~~~~~~~~~~~~ PREPARE OPERATOR ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = grady[exp(ik|x-y|)/|x-y|]
Gxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
Gxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
Gxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Finite elements
u = fem(sphere,typ);
v = fem(sphere,typ);

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,v);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy 
tic
Dbnd = 1/(4*pi) .* (integral(sigma,sigma,u,Gxy,ntimes(v),tol));
toc

% Regularization
tic
Dr   = 1/(4*pi) .* regularize(sigma,sigma,u,'grady[1/r]',ntimes(v));
Dbnd = Dbnd + Dr;
toc

% Operator [-Id/2 + Dt]
LHS = - 0.5*Id + Dbnd.';

% Finite element incident wave trace --> \int_Sx psi(x) dnx(pw(x)) dx
RHS = - integral(sigma,ntimes(u),gradxPW);


%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% LU factorization
tic
[Lh,Uh] = lu(LHS);
toc

% Solve linear system [-Id/2 + Dt] * lambda = dnP0
tic
lambda = Uh \ (Lh \ RHS); % LHS \ RHS;
toc


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
theta = 2*pi/1e3 .* (1:1e3)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) 1/(4*pi) .* exp(-1i*k*xdoty(X,Y));

% Finite element infinite operator --> \int_Sy exp(ik*nu.y) * psi(y) dx
Sinf = integral(nu,sigma,Ginf,v);

% Finite element radiation  
sol = Sinf * lambda;

% Analytical solution
ref = sphereHelmholtz('inf','neu',1,k,nu); 
norm(ref-sol,2)/norm(ref,2)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation
figure
plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')



%%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Green kernel
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dy 
tic
Sbnd = 1/(4*pi) .* integral(sigma,sigma,u,Gxy,v,tol);
toc

% Regularization
tic
Sbnd = Sbnd + 1/(4*pi) .* regularize(sigma,sigma,u,'[1/r]',v);
toc

% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy 
tic
Sdom = 1/(4*pi) .* integral(square.vtx,sigma,Gxy,v,tol);
toc

% Regularization
tic
Sdom = Sdom + 1/(4*pi) .* regularize(square.vtx,sigma,'[1/r]',v);
toc

% Boundary solution
Ibnd = integral(sigma,u,v);
Psca = Ibnd \ (Sbnd * lambda) ;
Pinc = PW(u.dof);
Pbnd = Pinc + Psca;

% Domain solution
Psca = Sdom * lambda;
Pinc = PW(square.vtx);
Pdom = Pinc + Psca;

% Annulation sphere interieure
r              = sqrt(sum(square.vtx.^2,2));
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
Pbnd = sphereHelmholtz('dom','neu',1,k,1.001*sphere.vtx) + PW(sphere.vtx);
Pdom = sphereHelmholtz('dom','neu',1,k,square.vtx) + PW(square.vtx);

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

