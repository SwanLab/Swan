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
%|    #    |   FILE       : nrtHmxHelmholtzBWneu.m                        |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Solve neumann scatering problem with          |
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
% camlight
% material dull
% lighting phong


%%% PREPARE OPERATORS
disp('~~~~~~~~~~~~~ PREPARE OPERATORS ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy         = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Finite elements
u = fem(sphere,typ);
v = fem(sphere,typ);

% Coupling coeff
beta = 1i*k;

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,v);

% Finite element boundary operator --> 
% k^2 * \int_Sx \int_Sy n.psi(x) G(x,y) n.psi(y) dx dy 
% - \int_Sx \int_Sy nxgrad(psi(x)) G(x,y) nxgrad(psi(y)) dx dy 
tic
Hr = 1/(4*pi) .* (k^2 * regularize(sigma,sigma,ntimes(u),'[1/r]',ntimes(v)) ...
    - regularize(sigma,sigma,nxgrad(u),'[1/r]',nxgrad(v)));
H  = 1/(4*pi) .* (k^2 * integral(sigma,sigma,ntimes(u),Gxy,ntimes(v),tol) ...
    - integral(sigma,sigma,nxgrad(u),Gxy,nxgrad(v),tol)) + Hr;
toc

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy 
tic
Dtr = 1/(4*pi) .* regularize(sigma,sigma,u,'grady[1/r]',ntimes(v)).';
Dt  = 1/(4*pi) .* integral(sigma,sigma,u,gradyGxy,ntimes(v),tol).' + Dtr;
toc

% Final operator [1i*k*beta*(-Id/2 + Dt) - H]
tic
LHS  = beta.*(-0.5*Id + Dt) - H;
toc

% Finite element incident wave trace --> \int_Sx psi(x) dnx(pw(x)) dx
RHS = - integral(sigma,ntimes(u),gradxPW);


%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% Preconditionneur ILU
tic
[L,U] = ilu(beta.*(-0.5*Id + Dtr) - Hr);
toc

% Solve linear system : [H + 1i*k*beta*(Id/2 - Dt)] = -dnP0
tic
mu = gmres(@(V) LHS*V,RHS,[],tol,100,L,U);
toc

% Jump for derivative
lambda = beta * mu;


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
theta = 2*pi/1e3 .* (1:1e3)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty        = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf         = @(X,Y) 1/(4*pi) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{1} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,1)) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{2} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,2)) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{3} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,3)) .* exp(-1i*k*xdoty(X,Y));

% Finite element infinite operators
Sinf = integral(nu,sigma,Ginf,v);
Dinf = integral(nu,sigma,gradxGinf,ntimes(v));

% Finite element radiation  
sol = Sinf*lambda - Dinf*mu;

% Analytical solution
ref = sphereHelmholtz('inf','neu',1,k,nu); 
norm(ref-sol,2)/norm(ref,2)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation
figure
plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')


%%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy         = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,v);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sbnd = 1/(4*pi) .* (integral(sigma,sigma,u,Gxy,v,tol) + ...
    regularize(sigma,sigma,u,'[1/r]',v));
toc

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy 
tic
Dbnd = 1/(4*pi) .* (integral(sigma,sigma,u,gradyGxy,ntimes(v),tol) + ...
    regularize(sigma,sigma,u,'grady[1/r]',ntimes(v)));
toc

% Boundary solution
Psca = Id\(Sbnd*lambda - (0.5*Id*mu + Dbnd*mu));
Pinc = PW(u.dof);
Pbnd = Pinc + Psca;

% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy 
tic
Sdom = 1/(4*pi) .* (integral(square.vtx,sigma,Gxy,v,tol) + ...
    regularize(square.vtx,sigma,'[1/r]',v));
toc

% Finite element radiative operator --> \int_Sx \int_Sy psi(x)' grady(G(x,y)) ny.psi(y) dx dy 
tic
Ddom = 1/(4*pi) .* ( integral(square.vtx,sigma,gradyGxy,ntimes(v),tol) + ...
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
