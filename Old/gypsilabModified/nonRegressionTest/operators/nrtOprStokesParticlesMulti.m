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
%|    #    |   FILE       : nrtOprStokesParticles.m                       |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Aline Lefevre-Lepot         |
%|  ( # )  |   CREATION   : 31.10.2018                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2018                                    |
%| ( === ) |   SYNOPSIS   : Movement of particles in a Stoke flow with    |
%|  `---'  |                obstacles                                     |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parametres
Nvtx = 100;
Nsph = 10;
tol  = 1e-3;
U0   = [1 0 0];
L    = 10;
Nprt = 100;
dt   = 0.01;  
Tf   = 20; 

% Random system of 3D particles
ctr   = -L/2+L*rand(Nsph,3);
[~,D] = knnsearch(ctr,ctr,'K',2);
rad   = (0.1+0.8*rand(Nsph,1)) .* (D(:,2)/2);
 
% Sphere mesh
sphere = mshSphere(Nvtx,1);

% Particles mesh
tic
vtx = cell(Nsph,1);
elt = cell(Nsph,1);
for i = 1:Nsph
    vtx{i} = ctr(i,:) + rad(i).*sphere.vtx;
    elt{i} = (i-1)*Nvtx + sphere.elt;
end
mesh = msh(cell2mat(vtx),cell2mat(elt));
toc

% Visu mesh
square        = mshSquare(1e3,[L L]);
visu          = square;
visu.vtx(:,3) = -L/2;
square.vtx    = [-L/2*ones(size(square.vtx,1),1) square.vtx(:,[1 2])];
visu          = union(visu,square);
square.vtx    = [square.vtx(:,3) L/2*ones(size(square.vtx,1),1) square.vtx(:,2)];
visu          = union(visu,square);

% Graphical representation
figure
plot(mesh)
hold on
plot(visu)
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
view(30,30)

% Finite element
gamma = dom(mesh,3); 
phi   = fem(mesh,'P1');

% Number of pool
Nlab = length(Composite);

% Domain decomposition for u
[Ilab,gammaLab,phiLab] = femSubdivide(gamma,phi,Nlab,0);
drawnow

% Stokeslet
tic
spmd
    M = cell(1,numlabs);
    for j = 1:Nlab
        % Operateur
        M{j} = oprIntegral('G',[],gammaLab{labindex},phiLab{labindex},...
            gammaLab{j},phiLab{j},tol);
        
        % Regularization
        M11  = regularize(gammaLab{labindex},gammaLab{j},phiLab{labindex},'[ij/r+rirj/r^3]11',phiLab{j});
        M12  = regularize(gammaLab{labindex},gammaLab{j},phiLab{labindex},'[ij/r+rirj/r^3]12',phiLab{j});
        M13  = regularize(gammaLab{labindex},gammaLab{j},phiLab{labindex},'[ij/r+rirj/r^3]13',phiLab{j});
        M22  = regularize(gammaLab{labindex},gammaLab{j},phiLab{labindex},'[ij/r+rirj/r^3]22',phiLab{j});
        M23  = regularize(gammaLab{labindex},gammaLab{j},phiLab{labindex},'[ij/r+rirj/r^3]23',phiLab{j});
        M33  = regularize(gammaLab{labindex},gammaLab{j},phiLab{labindex},'[ij/r+rirj/r^3]33',phiLab{j});
        M{j} = M{j} + 1/(8*pi) .* [M11 M12 M13 ; M12 M22 M23 ; M13 M23 M33];
    end
end
toc

% Convert to block matrix
tic
n = length(phi);
for i = 1:length(Ilab)
    Ilab{i} = [Ilab{i};n+Ilab{i};2*n+Ilab{i}];    
end
toc

% Define LHS
LHS = @(V) spmdProduct(Ilab,M,V);
        
% Boundary condition
muFct = {@(X) U0(1)*ones(size(X,1),1),...
    @(X) U0(2)*ones(size(X,1),1),...
    @(X) U0(3)*ones(size(X,1),1)};
mu = integral(gamma,phi,muFct);
mu = cell2mat(mu');

% Exact solver
% tic
% G = bmm(Ilab,Ilab,M);
% lambda = G\(-mu);
% toc
% figure
% spy(G)

% Iteractive solver
tic
lambda = mgcr(LHS,-mu,[],tol,100);
toc

% Reshape data
lambda = reshape(lambda,length(phi),3);

% Radiation
tic
U = zeros(size(visu.vtx,1),3);
for i=1:3
    for j = 1:3
        name   = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
        green  = @(X,Y) 1/(8*pi) .* femGreenKernel(X,Y,name,[]);
        U(:,i) = U(:,i) + integral(visu.vtx,gamma,green,phi,tol) * lambda(:,j);
    end
end
toc

% Add U0
U = U + U0;

% Graphical representation
for i = 1:3
    figure
    plot(mesh,'w')
    hold on
    plot(visu,U(:,i))
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    colorbar
    view(30,30)
end
drawnow


disp('~~> Michto gypsilab !')



