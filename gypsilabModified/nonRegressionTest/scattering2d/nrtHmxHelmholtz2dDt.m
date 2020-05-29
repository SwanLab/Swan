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
%|    #    |   FILE       : nrtHmxHelmholtz2dDt.m                         |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Martin Averseng             |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
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
X0  = [0 -1 0]
k   = 5

% Boundary mesh
mesh = mshCircle(N,1);

% Radiative mesh
radiat = mshSquare(10*N,[5 5]);

% Mesh representation
figure
plot(mesh)
hold on
plot(radiat)
plotNrm(mesh)
axis equal
axis(2.5*[-1 1 -1 1 -1 1])

% Frequency adjusted to maximum esge size
stp  = mesh.stp;
kmax = 1/stp(2)
% k    = kmax
if (k > kmax)
    warning('Wave number is too high for mesh resolution')
end
f = (k*340)/(2*pi);

% Incident wave
PW = @(X) exp(1i*k*X*X0');
gradxPW{1} = @(X) 1i*k*X0(1) .* PW(X);
gradxPW{2} = @(X) 1i*k*X0(2) .* PW(X);
gradxPW{3} = @(X) 1i*k*X0(3) .* PW(X);

% Incident wave representation
plot(radiat,real(PW(radiat.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
alpha(0.80)
view(0,90)


%%% PREPARE OPERATOR
disp('~~~~~~~~~~~~~ PREPARE OPERATOR ~~~~~~~~~~~~~')

% Green kernels
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]1',k);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]2',k);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]3',k);

% Domain
sigma = dom(mesh,gss);

% Finite elements
u = fem(mesh,typ);
v = fem(mesh,typ);

% Mass matrix
Id = integral(sigma,u,v);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dnyG(x,y) psi(y) dx dy 
tic
D = (1i/4) .* integral(sigma,sigma,u,dyGxy,ntimes(v),tol);
toc

% Regularization
tic
Dr = -1/(2*pi) .* regularize(sigma,sigma,u,'grady[log(r)]',ntimes(v));
toc

% Operator [-Id/2 + Dt]
LHS = - 0.5*Id + (D + Dr).';

% Finite element incident wave trace --> \int_Sx psi(x) dnx(pw(x)) dx
RHS = - integral(sigma,ntimes(u),gradxPW);

% Structure
figure
subplot(2,2,1:2)
spy(LHS)


%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% LU factorization
tic
[Lh,Uh] = lu(LHS);
toc
subplot(2,2,3)
spy(Lh)
subplot(2,2,4)
spy(Uh)

% Solve linear system [-Id/2 + Dt] * lambda = dnP0
tic
lambda = Uh \ (Lh \ RHS); % LHS \ RHS;
toc


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
theta = 2*pi/1e3 .* (1:1e3)';
nu    = [sin(theta),cos(theta),zeros(size(theta))];

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2); 
Ginf  = @(X,Y) exp(-1i*k*xdoty(X,Y));

% Finite element infinite operator --> \int_Sy exp(ik*nu.y) * psi(y) dx
Sinf = 1i/4 * integral(nu,sigma,Ginf,v);

% Finite element radiation  
sol = Sinf * lambda;

% Analytical solution
ref = diskHelmholtz('inf','neu',1,k,nu); 
norm(ref-sol,2)/norm(ref,2)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation
figure
plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')


%%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Green kernels
Gxy = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k);

% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy 
tic
Sdom = 1i/4 .* integral(radiat.vtx,sigma,Gxy,v,tol);
toc

% Regularization
tic
Sreg = -1/(2*pi) .* regularize(radiat.vtx,sigma,'[log(r)]',v);
Sdom = Sdom + Sreg;
toc

% Domain solution
Psca = Sdom * lambda;
Pinc = PW(radiat.vtx);
Pdom = Psca + Pinc;

% Annulation mesh interieure
r             = sqrt(sum(radiat.vtx.^2,2));
Pdom(r<=1.01) = Pinc(r<=1.01);

% Graphical representation
figure
plot(radiat,abs(Pdom))
axis equal
title('Total field solution')
colorbar


%%% ANAYTICAL SOLUTIONS FOR COMPARISONS
% Analytical solution
Pdom = diskHelmholtz('dom','neu',1,k,radiat.vtx) + PW(radiat.vtx);

% Solution representation
figure
plot(radiat,abs(Pdom))
axis equal;
title('Analytical solution')
colorbar
view(0,90)



disp('~~> Michto gypsilab !')


