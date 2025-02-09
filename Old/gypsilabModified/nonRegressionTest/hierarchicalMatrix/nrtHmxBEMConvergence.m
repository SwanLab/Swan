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
%|    #    |   FILE       : nrtHmxBEMConvergence.m                        |
%|    #    |   VERSION    : 0.60                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Solve dirichlet scatering problem with single |
%|  `---'  |                layer potential using H-Matrix solver         |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N2  = 2.^(10:13)
tol = 1e-6
typ = 'P1'
gss = 3
X0  = [0 0 -1]
k   = 5
f   = (k*340)/(2*pi)

% Incident wave
PW = @(X) exp(1i*k*X*X0');

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Plane waves direction
theta = 2*pi/1e3 .* (1:1e3)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Analytical solution
ref = sphereHelmholtz('inf','dir',1,k,nu);

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) 1/(4*pi) .* exp(-1i*k*xdoty(X,Y));

% Loop for N
timeBld = zeros(size(N2));
timeSol = zeros(size(N2));
err     = zeros(size(N2));
for i = 1:length(N2)
    % Mesh size
    N = N2(i)
    
    % Spherical mesh
    sphere = mshSphere(N,1);
    
    % Domain
    sigma = dom(sphere,gss);
    
    % Graphical representation
    figure(10); clf
    plot(sphere)
    axis equal
    plot(sphere,real(PW(sphere.vtx)))
    title('Incident wave')
    xlabel('X');   ylabel('Y');   zlabel('Z');
    hold off
    view(0,10)
    
    % Finite elements
    u = fem(sphere,typ);
    v = fem(sphere,typ);
    
    % Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy
    tic
    LHS = 1/(4*pi) .* integral(sigma,sigma,u,Gxy,v,tol);
    
    % Regularization
    Sr  = 1/(4*pi) .* regularize(sigma,sigma,u,'[1/r]',v);
    LHS = LHS + Sr;
    
    % Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
    RHS = - integral(sigma,u,PW);
    timeBld(i) = toc
    
    % Solve linear system [S] * lambda = P0
    tic
    lambda     = LHS \ RHS;
    timeSol(i) = toc
    
    % Finite element infinite operator --> \int_Sy exp(ik*nu.y) * psi(y) dx
    Sinf = integral(nu,sigma,Ginf,v);
    
    % Finite element radiation
    sol = Sinf * lambda;
    
    % Error
    % err = norm(ref-sol,2)/norm(ref,2)
    err(i) = norm(ref-sol,'inf')/norm(ref,'inf')
    
    % Graphical representation
    figure(11); clf
    plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')
end

% Graphical representation
figure
loglog(N2,err,'b+-')
hold on
loglog([1 10]*1e3,[1 1/10]*1e-3,'r')
grid on
xlabel('Nddl')
ylabel('max(|Pref-Psol|)/max(|Pref|)')
legend({'Error','Slope -1'})
title('Relative error for infinite field')

figure
loglog(N2,timeBld,'b-+')
hold on
loglog([1 10]*1e3,[1 10],'r')
grid on
xlabel('Nddl')
ylabel('Time (s)')
title('Assembling time (H-Matrix 10^{-6})')

figure
loglog(N2,timeSol,'b-+')
hold on
loglog([1 10]*1e3,[1 10],'r')
grid on
xlabel('Nddl')
ylabel('Time (s)')
title('Exact solver time  (H-Matrix 10^{-6})')




disp('~~> Michto gypsilab !')


