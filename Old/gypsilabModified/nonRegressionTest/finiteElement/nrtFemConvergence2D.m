%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Francois Alouges (c) 2017-2018.                            |
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
%|    #    |   FILE       : nrtFemConvergence2D.m                         |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Francois Alouges                              |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Mesh convergence for lagrange approximation   |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

errL2P1 = []; errH1P1 = errL2P1; errL2P2 = []; errH1P2 = errL2P2; h = [];
Uex = @(x) cos(x(:,1)).*cos(x(:,2))/3;
f = @(x) cos(x(:,1)).*cos(x(:,2));

for nbPts=10:3:80
    % Square mesh
    mesh  = mshSquare(nbPts^2,[2*pi 2*pi]);
    omega = dom(mesh,7);
    h     = [h, 2*pi/nbPts];
%     plot(omega);
    
    
    %%% SOLVE LINEAR PROBLEM
    disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')
    
    
    % Finite elements
    u  = fem(mesh,'P1'); v  = fem(mesh,'P1');
    u2 = fem(mesh,'P2'); v2 = fem(mesh,'P2');
    
    % Matrix and RHS
    K  = integral(omega,grad(u),grad(v)) + integral(omega,u,v);
    K2 = integral(omega,grad(u2),grad(v2)) + integral(omega,u2,v2);
    F  = integral(omega, v, f);
    F2 = integral(omega, v2, f);
    
    % Résolution
    uh = K\F;
    uh2 = K2\F2;

    % erreur en norme L2 et H1
    errL2P1 = [errL2P1, omega.diff(u, uh, Uex, 'L2')];
    errH1P1 = [errH1P1, omega.diff(u, uh, Uex,  'H1')];
    errL2P2 = [errL2P2, omega.diff(u2, uh2, Uex,  'L2')];
    errH1P2 = [errH1P2, omega.diff(u2, uh2, Uex,  'H1')];
end

% Plot the error graphs
figure(2)
subplot(1,2,1)
loglog(h,errL2P1,'b+-',h,errL2P2,'r+-',h,(h*10^-0.5).^2,'k--',h,(h*10^-0.5).^3,'k:')
legend('EF P1', 'EF P2','slope 2','slope 3')
title('error L2, CB Neumann')

xlabel('h')
subplot(1,2,2)
loglog(h,errH1P1,'b+-',h,errH1P2,'r+-',h,h,'k',h,(h*10^-0.5).^2,'k--')
legend('EF P1', 'EF P2','slope 1','slope 2')
title('error H1, CB Neumann')
xlabel('h')


disp('~~> Michto gypsilab !')

