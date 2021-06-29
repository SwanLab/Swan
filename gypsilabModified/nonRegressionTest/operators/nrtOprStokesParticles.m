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
tol  = 1e-3;
U0   = [1 0 0];
L    = 10;
Nprt = 100;
dt   = 0.01;  
Tf   = 20; 
 
% Sphere mesh
mesh = mshSphere(Nvtx,1);

% Visu mesh
square     = mshSquare(L*Nvtx,[L L]);
visu       = square;
square.vtx = square.vtx(:,[3 1 2]);
visu       = union(visu,square);
square.vtx = square.vtx(:,[2 1 3]);
visu       = union(visu,square);

% Graphical representation
figure
plot(mesh)
hold on
plot(visu)
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')

% Finite element
gamma = dom(mesh,3); 
phi   = fem(mesh,'P1');

% Number of pool
Nlab = length(Composite);

% Domain decomposition for u
[Ilab,gammaLab,phiLab] = femSubdivide(gamma,phi,Nlab,10);
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
G = bmm(Ilab,Ilab,M);
toc
        
% Graphical representation
figure
spy(G)

% Boundary condition
muFct = {@(X) U0(1)*ones(size(X,1),1),...
    @(X) U0(2)*ones(size(X,1),1),...
    @(X) U0(3)*ones(size(X,1),1)};
mu = integral(gamma,phi,muFct);
mu = cell2mat(mu');

% Resolution
tic
lambda = G\(-mu);
lambda = reshape(lambda,length(phi),3);
toc

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
    plot(mesh)
    hold on
    plot(visu,U(:,i))
    axis equal
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    colorbar
end
drawnow

% Particles coordinates
Xprt = [-L/2*ones(Nprt,1),-L/2+L*rand(Nprt,1),zeros(Nprt,1)];

% Particles representation
figure
plot(mesh)
hold on
plot(visu,sqrt(sum(U.^2,2)))
H = plot3(Xprt(:,1),Xprt(:,2),Xprt(:,3),'*w');
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
colorbar

% Boucle temporellle
t = 0;
while (t<Tf)
    % Plus proche voisin
    I = knnsearch(visu.vtx,Xprt,'K',1);
    
    % Deplacement
    Xprt = Xprt + dt*U(I,:);
    
    % Visu
    set(H,'XData',Xprt(:,1),'YData',Xprt(:,2),'ZData',Xprt(:,3))
    drawnow
    
    % Incrementation
    t = t + dt; 
end


disp('~~> Michto gypsilab !')



