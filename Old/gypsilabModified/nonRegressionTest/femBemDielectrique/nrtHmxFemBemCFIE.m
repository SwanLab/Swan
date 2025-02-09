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
%|    #    |   FILE       : nrtHmxFemBemCFIE.m                            |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Francois Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Solve fem/bem with PEC scatering problem      |
%|  `---'  |                Volumic ans surfacic CFIE formulation         |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')


%%% PREPARATION
disp('~~~~~~~~~~~~~ PREPARATION ~~~~~~~~~~~~~')

% Spherical mesh
mesh = msh('sphere_1e3.msh');

% Normalisation sphere unite interieure
mesh.vtx = mesh.vtx/0.8;

% Frontieres surfaciques
bound = mesh.bnd;
Iint  = (sum(bound.ctr.*bound.nrm,2) < 0); 
int   = bound.sub(Iint);
ext   = bound.sub(~Iint);

% Domaine volumique
omega = dom(mesh,4);

% Domaine exterieur
sigma = dom(ext,3);

% Graphical representation
figure
plot(int,'r')
hold on
plot(ext,'b')
alpha(0.5)
axis equal

% Frequency adjusted to maximum esge size
stp = mesh.stp;
k   = 1/stp(2);
c   = 299792458;
f   = (k*c)/(2*pi);
disp(['Frequency : ',num2str(f/1e6),' MHz']);

% Accuracy
tol = 1e-3;
disp(['Accuracy : ',num2str(tol)]);

% Coupling factor
beta = 0.5;
disp(['CFIE factor : ',num2str(beta)]);

% Incident direction and field
X0 = [0 0 -1]; 
E  = [0 1  0]; % Polarization (+x for Theta-Theta and +y for Phi-Phi)
H  = cross(X0,E);

% Incident Plane wave (electromagnetic field)
PWE{1} = @(X) exp(1i*k*X*X0') * E(1);
PWE{2} = @(X) exp(1i*k*X*X0') * E(2);
PWE{3} = @(X) exp(1i*k*X*X0') * E(3);

% Incident Plane wave (electromagnetic field)
PWH{1} = @(X) exp(1i*k*X*X0') * H(1);
PWH{2} = @(X) exp(1i*k*X*X0') * H(2);
PWH{3} = @(X) exp(1i*k*X*X0') * H(3);

% Incident wave representation
figure
plot(mesh,real(PWE{2}(mesh.vtx)))
alpha(0.1)
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
axis equal
view(0,10)


%%% FORMUMATIONS
disp('~~~~~~~~~~~~~ FORMULATIONS ~~~~~~~~~~~~~')

% Green kernel
Gxy    = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
Hxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k) ;
Hxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k) ;
Hxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k) ;
Ixy{1} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]1',k) ;
Ixy{2} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]2',k) ;
Ixy{3} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]3',k) ;

% Surfacic finite elements for magnetic current
Jh = fem(ext,'RWG');
N1 = size(Jh.unk,1);

% Volumique finite elements for electric field
Eh = fem(mesh,'NED');
Eh = dirichlet(Eh,int);
N2 = size(Eh.unk,1);

% Restriction matrix for regularization
[~,I1,I2] = intersect(Jh.unk,Eh.unk,'rows');
P         = sparse(I1,I2,1,size(Jh.unk,1),size(Eh.unk,1));

% Vectors
JEi   = integral(sigma,Jh,PWE);
JnxHi = integral(sigma,nx(Jh),PWH);

% Sparse operators
tic
rotErotE = integral(omega,curl(Eh),curl(Eh));
EE       = integral(omega,Eh,Eh);
JE       = integral(sigma,Jh,Eh);
JJ       = integral(sigma,Jh,Jh);
toc

% Full operators
tic
JTJ  = 1/(4*pi) .* integral(sigma, sigma, Jh, Gxy, Jh,tol) ...
    - 1/(4*pi*k^2) * integral(sigma, sigma, div(Jh), Gxy, div(Jh),tol) ;
JTJr = 1/(4*pi) * regularize(sigma, sigma, Jh, '[1/r]', Jh) ...
      - 1/(4*pi*k^2) * regularize(sigma, sigma, div(Jh), '[1/r]', div(Jh));
JTJ = JTJ + JTJr;
toc

tic
JnxKJ  = 1/(4*pi) * integral(sigma, sigma, nx(Jh), Hxy, Jh, tol); 
JnxKJr = 1/(4*pi) * regularize(sigma, sigma, nx(Jh), 'grady[1/r]', Jh);
JnxKJ  = JnxKJ + JnxKJr;
toc

tic
JKExn  = - 1/(4*pi) * integral(sigma, sigma, Jh, Hxy, nx(Eh),tol); 
JKExnr = - 1/(4*pi) * regularize(sigma, sigma, Jh, 'grady[1/r]', Jh) * P; 
JKExn  = JKExn + JKExnr;
toc

tic
JnxTExn  = 1/(4*pi) .* integral(sigma, sigma, nx(Jh), Gxy, nx(Eh),tol) ...
    + 1/(4*pi*k^2) * integral(sigma, sigma, nx(Jh), Ixy, divnx(Eh),tol);
JnxTExnr = 1/(4*pi) * regularize(sigma, sigma, nx(Jh), '[1/r]', Jh) * P ;
JnxTExn  = JnxTExn + JnxTExnr;
toc

% LHS = [A B ; C D]
A  = beta * (1i*k*JTJ) + (1-beta) * (0.5*JJ - JnxKJ);
Ar = beta * (1i*k*JTJr) + (1-beta) * (0.5*JJ - JnxKJr);
size(A)

B  = beta * (JKExn - 0.5*JE) + (1-beta) * (-1i*k*JnxTExn);
Br = beta * (JKExnr - 0.5*JE) + (1-beta) * (-1i*k*JnxTExnr);
size(B)

C = JE.';
size(C)

D = 1/(1i*k) * (rotErotE - k^2*EE);
size(D)

% Right hand side
Y   = beta * (-JEi) + (1-beta) * (-JnxHi); 
RHS = [Y;zeros(size(D,1),1)];


%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% Factorization LU H-Matrix
tic
[L,U] = ilu(Ar);
toc

% Shurr complement resolution
tic
SV   = @(V) A*V - B*(D\(C*V));
J    = gmres(SV,Y,[],tol,100,L,U);
E    = - D\(C*J);
toc


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
Ninf  = 1e3;
theta = 2*pi/1e3 .* (1:Ninf)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) exp(-1i*k*xdoty(X,Y));

% Magnetic current radiation
Tinf = integral(nu,sigma,Ginf,Jh);
Jinf = 1i*k/(4*pi)*cross(nu, cross([Tinf{1}*J, Tinf{2}*J, Tinf{3}*J], nu));

% Electric field radiation 
Kinf = integral(nu,sigma,Ginf,nx(Eh));
Einf = 1i*k/(4*pi) * cross(nu, [-Kinf{1}*E, -Kinf{2}*E, -Kinf{3}*E] );

% Total electric field radiation
sol = Jinf - Einf;   

% Radiation infinie de reference, convention e^(+ikr)/r
nMax = 100; refInf = zeros(Ninf,1);
if E(1) == 1
    for jj = 1:Ninf
        refInf(jj,:) = sphereMaxwell(1, -f, theta(jj), 0.0, nMax);
    end
else
    for jj = 1:Ninf
        [~,refInf(jj)] = sphereMaxwell(1, -f, theta(jj), pi/2, nMax);
    end
end
refInf = refInf ./ sqrt(4*pi);

% Radiations infinies en X
if E(1) == 1
    sol = sin(theta)'.*sol(:,3) - cos(theta)'.*sol(:,1);
else
    sol = sol(:,2);
end

% Erreur
eL2   = norm(refInf-sol,2)/norm(refInf,2)
eLINF = norm(refInf-sol,'inf')/norm(refInf,'inf')
             
% Representation graphique
figure
subplot(1,2,1)
plot(theta,20*log10(abs(sol)),'b',theta,20*log10(abs(refInf)),'r--')

subplot(1,2,2)
plot(theta,real(sol),'--b', theta,imag(sol),'--r', theta, real(refInf),':b', theta,imag(refInf),':r');
drawnow



disp('~~> Michto gypsilab !')


