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
%|    #    |   FILE       : nrtOprSpmd.m                                  |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2018                                    |
%| ( === ) |   SYNOPSIS   : Parallel build for operators                  |
%|  `---'  |                                                              |
%+========================================================================+

% Nettoyage
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Dimensions
N    = 5e2;
tol  = 1e-3;
k    = 5;
Nlab = length(Composite);

% Meshes
mesh = mshSphere(N,1);

% Domain
sigma = dom(mesh,3);

% Finite element
u = fem(mesh,'P1');

% Particles charges
V = (-1+2*rand(length(u),1)) + (-1+2i*rand(length(u),1));

% Spatial representation of particles
figure
plot(mesh,'b')
alpha(0.5)
axis equal 

% Domain decomposition for u
[Ilab,sigmaLab,uLab] = femSubdivide(sigma,u,Nlab,10);
drawnow

% Parallel
tic
spmd
    M = cell(1,numlabs);
    for j = 1:Nlab   
        M{j} = oprIntegral('G',k,sigmaLab{labindex},uLab{labindex},...
            sigmaLab{j},uLab{j},tol);
    end
end
toc

% Sequential
tic
Mref = oprIntegral('G',k,sigma,u,sigma,u,tol);
toc

% Special issue Stokes
n = length(u);
for i = 1:length(Ilab)
    Ilab{i} = [Ilab{i};n+Ilab{i};2*n+Ilab{i}];    
end
V = (-1+2*rand(3*length(u),2)) + (-1+2i*rand(3*length(u),2));

% Convert to block matrix
Msol = bmm(Ilab,Ilab,M);

% Matrix-vector product
tic
ref = Mref * V;
toc
tic
sol = Msol * V;
toc
tic
sol2 = spmdProduct(Ilab,M,V);
toc

% Error
norm(ref-sol,'inf')/norm(ref,'inf')
norm(ref-sol2,'inf')/norm(ref,'inf')

% Graphical representation
figure
spy(Mref)
figure
spy(Msol)




disp('~~> Michto gypsilab !')


