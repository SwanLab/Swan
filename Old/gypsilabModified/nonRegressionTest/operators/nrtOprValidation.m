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
%|    #    |   FILE       : nrtOprValidation.m                            |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Operators validation                          |
%|  `---'  |                                                              |
%+========================================================================+

% Nettoyage
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensions
Nx = 1e3;
Ny = 5e2;

% Accuracy
tol = 1e-3;

% Meshes
Xmsh = mshSphere(Nx,1);
Ymsh = mshSphere(Ny,1);%Xmsh;%mshCube(Ny,2*[1 1 1]);
% Ymsh = Ymsh.bnd;

% Domain
Xdom = dom(Xmsh,3);
Ydom = dom(Ymsh,3);

% Wave number or frequency (Hz)
k = 5;

% Green kernels 
Gxy         = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
gradxGxy{1} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]1',k);
gradxGxy{2} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]2',k);
gradxGxy{3} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]3',k);
gradyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
gradyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
gradyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Spatial representation of particles
figure
plot(Xmsh,'b')
hold on
plot(Ymsh,'r')
alpha(0.5)
axis equal 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELMHOLTZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ HELMHOLTZ ~~~~~~~~~~~~~')

% Finite element
Xfem = fem(Xmsh,'P1');
Yfem = fem(Ymsh,'P1');

% Particles charges
V = (-1+2*rand(length(Yfem),2)) + (-1+2i*rand(length(Yfem),2));

% Single layer
tic
Mref = 1/(4*pi) .* integral(Xdom,Ydom,Xfem,Gxy,Yfem,tol);
toc
tic
Msol = oprIntegral('S',k,Xdom,Xfem,Ydom,Yfem,tol);
toc
norm(Mref*V-Msol*V,'inf')/norm(Mref*V,'inf')

% Double layer
tic
Mref = 1/(4*pi) .* integral(Xdom,Ydom,Xfem,gradyGxy,ntimes(Yfem),tol);
toc
tic
Msol = oprIntegral('D',k,Xdom,Xfem,Ydom,Yfem,tol);
toc
norm(Mref*V-Msol*V,'inf')/norm(Mref*V,'inf')

% Double layer transpose
tic
Mref = 1/(4*pi) .* integral(Xdom,Ydom,ntimes(Xfem),gradxGxy,Yfem,tol);
toc
tic
Msol = oprIntegral('Dt',k,Xdom,Xfem,Ydom,Yfem,tol);
toc
norm(Mref*V-Msol*V,'inf')/norm(Mref*V,'inf')

% Hypersingular
tic
Mref = 1/(4*pi) .* (k^2 .* integral(Xdom,Ydom,ntimes(Xfem),Gxy,ntimes(Yfem),tol) + ...
    - integral(Xdom,Ydom,nxgrad(Xfem),Gxy,nxgrad(Yfem),tol));
toc
tic
Msol = oprIntegral('H',k,Xdom,Xfem,Ydom,Yfem,tol);
toc
norm(Mref*V-Msol*V,'inf')/norm(Mref*V,'inf')

disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAXWELL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ MAXWELL ~~~~~~~~~~~~~')
tic

% Finite element
Xfem = fem(Xmsh,'RWG');
Yfem = fem(Ymsh,'RWG');

% Particles charges
V = (-1+2*rand(length(Yfem),2)) + (-1+2i*rand(length(Yfem),2));

% EFIE
tic
Mref = 1/(4*pi) .* (1i*k .* integral(Xdom,Ydom,Xfem,Gxy,Yfem,tol) + ...
    -1i/k .* integral(Xdom,Ydom,div(Xfem),Gxy,div(Yfem),tol));
toc
tic
Msol = oprIntegral('T',k,Xdom,Xfem,Ydom,Yfem,tol);
toc
norm(Mref*V-Msol*V,'inf')/norm(Mref*V,'inf')

% MFIE
tic
Mref = 1/(4*pi) * integral(Xdom,Ydom,nx(Xfem),gradyGxy,Yfem,tol);
toc
tic
Msol = oprIntegral('nxK',k,Xdom,Xfem,Ydom,Yfem,tol);
toc
norm(Mref*V-Msol*V,'inf')/norm(Mref*V,'inf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ STOKES ~~~~~~~~~~~~~')
tic

% Finite element
Xfem = fem(Xmsh,'P1');
Yfem = fem(Ymsh,'P1');

% Particles charges
V = (-1+2*rand(3*length(Yfem),2)) + (-1+2i*rand(3*length(Yfem),2));

% Symetry indices
ind = [1 1 ; 1 2 ; 1 3 ; 2 2 ; 2 3 ; 3 3];
M   = cell(6,1);
un  = ones(length(Yfem),1);

% Stokeslet
tic
for l = 1:length(ind)
    % Local indices
    i = ind(l,1);
    j = ind(l,2);
    
    % Green kernel
    name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
    green = @(X,Y) 1/(8*pi) .* femGreenKernel(X,Y,name,[]);
    M{l}  = integral(Xdom,Ydom,Xfem,green,Yfem,tol);
end
Mref = bmm({M{1}    M{2}    M{3}    ; M{2}    M{4}    M{5}    ; M{3}    M{5}    M{6}});
toc
tic
Msol = oprIntegral('G',k,Xdom,Xfem,Ydom,Yfem,tol);
toc
norm(Mref*V-Msol*V,'inf')/norm(Mref*V,'inf')

% Stresslet
tic
for l = 1:length(ind)
    % Local indices
    i = ind(l,1);
    j = ind(l,2);

    % Gren kernel
    M{l} = zeros(M{l});
%     P    = 0;
    for k = 1:3
        % Stresslet
        name  = ['[rirjrk/r^5]',num2str(i),num2str(j),num2str(k)];
        green = @(X,Y) -6/(8*pi) .* femGreenKernel(X,Y,name,[]);
        M{l}  = M{l} + integral(Xdom,Ydom,Xfem,green,ntimes(Yfem,k),tol);
        
        % Correction
%         P = P + integral(Xdom.qud,Ydom,green,ntimes(Yfem,k),tol) * un;
    end
    
    % Correction
%     M{l} = M{l} - integral(Xdom,Xfem,@(X)P,Yfem);
end
Mref = bmm({M{1}    M{2}    M{3}    ; M{2}    M{4}    M{5}    ; M{3}    M{5}    M{6}});
toc
tic
Msol = oprIntegral('Ts',k,Xdom,Xfem,Ydom,Yfem,tol);
toc
norm(Mref*V-Msol*V,'inf')/norm(Mref*V,'inf')





disp('~~> Michto gypsilab !')


