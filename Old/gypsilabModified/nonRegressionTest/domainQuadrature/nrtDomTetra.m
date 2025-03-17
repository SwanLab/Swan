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
%|    #    |   FILE       : nrtDomTetra.m                                 |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Tetrahedral quadrature and integration        |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nvtx = 1e2;
L    = [4 3 2];
tol  = 1e-3;

% Square mesh
mesh  = mshCube(Nvtx,L);
bound = mesh.bnd;

% Domain
omega = dom(mesh,4);

% Finite element space
u = fem(mesh,'P1');
v = fem(mesh,'P1');

% Graphical representation
figure
plot(mesh)
hold on
plot(omega)
plot(u,'or')
axis equal
alpha(0.1)
view(30,30)

% Numerical functions
Fx   = @(X) ones(size(X,1),1);
Fx3  = {Fx,Fx,Fx};
Fxy  = @(X,Y) (X(:,1)==Y(:,1)) .* (X(:,2)==Y(:,2)) .* (X(:,3)==Y(:,3)) ;
Fxy3 = {Fxy,Fxy,Fxy};
Gxy  = @(X,Y) 1./(4*pi) * femGreenKernel(X,Y,'[exp(ikr)/r]',5);
Gxy3 = {Gxy,Gxy,Gxy};



%%%%%%%%%%%%%%% 2 ARGUMENTS %%%%%%%%%%%%%%%
disp('=========== 2 ARGUMENTS  ============')

% \int_{mesh(x)} f(x) dx 
ref = integral(omega,Fx);
abs(ref-prod(L))/prod(L)

% \int_{mesh(x)} f3(x) dx 
sol = integral(omega,Fx3);
norm(sol{1}-ref,'inf')



%%%%%%%%%%%%%%% 3 ARGUMENTS %%%%%%%%%%%%%%%
disp('=========== 3 ARGUMENTS  ============')

% \int_{mesh(y)} f(x,y) dy 
ref = integral(omega.qud,omega,Fxy);
abs(sum(ref,1)-prod(L))/prod(L)

% \int_{mesh(y)} f3(x,y) dy 
sol = integral(omega.qud,omega,Fxy3);
norm(sol{1}-ref,'inf')

%-------------------------------------

% \int_{mesh(x)} f(x,y) dx 
ref = integral(omega,omega.qud,Fxy);
abs(sum(ref,2)-prod(L))/prod(L)

% \int_{mesh(y)} f3(x,y) dx 
sol = integral(omega,omega.qud,Fxy3);
norm(sol{1}-ref,'inf')

%-------------------------------------

% \int_{mesh(x)} f(x) psi(x) dx 
ref = integral(omega,Fx,v);
abs(sum(ref,2) - prod(L))/prod(L)

% \int_{mesh(x)} f3(x) psi(x) dx 
sol = integral(omega,Fx3,v);
norm(sol{1}-ref,'inf')

%-------------------------------------

% \int_{mesh(x)} psi(x)' f(x)  dx 
ref = integral(omega,u,Fx);
abs(sum(ref,1) - prod(L))/prod(L)

% \int_{mesh(x)} psi(x)' f3(x)  dx 
sol = integral(omega,u,Fx3);
norm(sol{1}-ref,'inf')

%-------------------------------------

% \int_{mesh(x)} psi(x)' psi(x) dx 
ref = integral(omega,u,v);
abs(sum(sum(ref,1),2) - prod(L))/prod(L)



%%%%%%%%%%%%%%% 4 ARGUMENTS %%%%%%%%%%%%%%%
disp('=========== 4 ARGUMENTS  ============')

% \int_{mesh(x)} psi(x)' f(x) psi(x) dx 
ref = integral(omega,u,Fx,v);
abs(sum(sum(ref,1),2) - prod(L))/prod(L)

% \int_{mesh(x)} psi(x)' f3(x) psi(x) dx 
sol = integral(omega,u,Fx3,v);
norm(sol{3}-ref,'inf')

%-------------------------------------

% \int_{mesh(y)} f(x,y) psi(y) dy    
ref = integral(omega.qud,omega,Fxy,v);
abs(sum(sum(ref,1),2) - prod(L))/prod(L)

% \int_{mesh(y)} f(3x,y) psi(y) dy    
sol = integral(omega.qud,omega,Fxy3,v);
norm(sol{3}-ref,'inf')

%-------------------------------------

% \int_{mesh(x)} psi(x)' f(x,y) dx  
ref = integral(omega,omega.qud,u,Fxy);
abs(sum(sum(ref,1),2) - prod(L))/prod(L)

% \int_{mesh(x)} psi(x)' f3(x,y) dx  
sol = integral(omega,omega.qud,u,Fxy3);
norm(sol{3}-ref,'inf')



%%%%%%%%%%%%%% 5 ARGUMENTS %%%%%%%%%%%%%%%
disp('=========== 5 ARGUMENTS  ============')

% \int_{mesh(y)} G(x,y) psi(y) dy    
ref = integral(mesh.ctr,omega,Gxy,v);
sol = integral(mesh.ctr,omega,Gxy,v,tol);
norm(full(sol)-ref,'inf')./norm(ref,'inf')

% % \int_{mesh(y)} G3(x,y) psi(y) dy    
sol = integral(mesh.ctr,omega,Gxy3,v,tol);
norm(full(sol{1})-ref,'inf')./norm(ref,'inf')

%-------------------------------------

% \int_{mesh(x)} psi(x)' G(x,y) dx    
ref = integral(omega,mesh.ctr,u,Gxy);
sol = integral(omega,mesh.ctr,u,Gxy,tol);
norm(full(sol)-ref,'inf')./norm(ref,'inf')

% \int_{mesh(x)} psi(x)' G3(x,y) dx    
sol = integral(omega,mesh.ctr,u,Gxy3,tol);
norm(full(sol{1})-ref,'inf')./norm(ref,'inf')



%%%%%%%%%%%%%%% 6 ARGUMENTS %%%%%%%%%%%%%%%
disp('=========== 6 ARGUMENTS  ============')

% \int_{mesh(x)} \int_{mesh(y)} psi(x)' G(x,y) psi(y) dx dy    
ref = integral(omega,omega,u,Gxy,v);
sol = integral(omega,omega,u,Gxy,v,tol);
norm(full(sol)-ref,'inf')./norm(ref,'inf')

% \int_{mesh(x)} \int_{mesh(y)} grad(psi(x))' G(x,y) psi(y) dx dy    
ref = integral(omega,omega,grad(u),Gxy,v);
sol = integral(omega,omega,grad(u),Gxy,v,tol);
norm(full(sol{3})-ref{3},'inf')./norm(ref{3},'inf')

% \int_{mesh(x)} \int_{mesh(y)} psi(x)' G3(x,y) psi(y) dx dy    
ref = integral(omega,omega,u,Gxy3,v);
sol = integral(omega,omega,u,Gxy3,v,tol);
norm(full(sol{3})-ref{3},'inf')./norm(ref{3},'inf')

% \int_{mesh(x)} \int_{mesh(y)} psi(x)' G(x,y) grad(psi(y)) dx dy    
ref = integral(omega,omega,u,Gxy,grad(v));
sol = integral(omega,omega,u,Gxy,grad(v),tol);
norm(full(sol{3})-ref{3},'inf')./norm(ref{3},'inf')

% \int_{mesh(x)} \int_{mesh(y)} grad(psi(x))' G3(x,y) psi(y) dx dy    
ref = integral(omega,omega,grad(u),Gxy3,v);
sol = integral(omega,omega,grad(u),Gxy3,v,tol);
norm(full(sol)-ref,'inf')./norm(ref,'inf')

% \int_{mesh(x)} \int_{mesh(y)} psi(x)' G3(x,y) grad(psi(y)) dx dy    
ref = integral(omega,omega,u,Gxy3,grad(v));
sol = integral(omega,omega,u,Gxy3,grad(v),tol);
norm(full(sol)-ref,'inf')./norm(ref,'inf')

% \int_{mesh(x)} \int_{mesh(y)} grad(psi(x))' G(x,y) grad(psi(y)) dx dy    
ref = integral(omega,omega,grad(u),Gxy,grad(v));
sol = integral(omega,omega,grad(u),Gxy,grad(v),tol);
norm(full(sol)-ref,'inf')./norm(ref,'inf')

% \int_{mesh(x)} \int_{mesh(y)} grad(psi(x))' G3(x,y) grad(psi(y)) dx dy    
ref = integral(omega,omega,grad(u),Gxy3,grad(v));
sol = integral(omega,omega,grad(u),Gxy3,grad(v),tol);
norm(full(sol)-ref,'inf')



disp('~~> Michto gypsilab !')




