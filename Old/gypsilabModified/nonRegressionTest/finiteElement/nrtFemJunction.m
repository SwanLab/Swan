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
%|    #    |   FILE       : nrtFemDirichletSquare.m                       |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Dirichlet condition with a square             |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nvtx = 1e3;
Neig = 10;
L    = [1 0.5];
dl   = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREPARE MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh1          = mshSquare(ceil(Nvtx/2),[L(1)/2 L(2)]);
% mesh1.vtx(:,1) = mesh1.vtx(:,1) - L(1)/4;
% mesh1.col(:)    = 1;
% 
% mesh2          = mesh1;
% mesh2.vtx(:,1) = mesh2.vtx(:,1) + L(1)/2;
% mesh2.col(:)   = 2;
% 
% mesh3          = mesh2;
% mesh3.vtx(:,3) = mesh3.vtx(:,1);
% mesh3.vtx(:,1) = 0;
% ctr            = mesh3.ctr;
% mesh3          = mesh3.sub(ctr(:,2) <= 0.13);
% ctr            = mesh3.ctr;
% mesh3          = mesh3.sub(ctr(:,2) >= -0.13);
% ctr            = mesh3.ctr;
% mesh3          = mesh3.sub(ctr(:,3) <= 0.07);
% mesh3.col(:)   = 3;
% 
% mesh  = union(union(mesh1,mesh2),mesh3);
% mshWriteVtk('Tjunction.vtk',mesh,mesh.col)
% pouet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read mesh and boundary
mesh  = msh('junction.vtk');
bound = mesh.bnd;

% Extract Tjunction
ctr       = mesh.ctr;
Tjunction = mesh.sub(ctr(:,3)>0);

% Extract interface
int = setdiff(Tjunction.bnd,bound);

% Extract vextex cut
cut = int.sub(2:length(int)-1);

% Extract plaque 
plaque = setdiff(mesh,Tjunction);

% Subdivide plaque
ctr     = plaque.ctr;
plaque1 = plaque.sub(ctr(:,1)<=0);
plaque2 = plaque.sub(ctr(:,1)>=0);

% Cut and translate plaque 1
I = ismember(plaque1.vtx,cut.vtx,'rows');
plaque1.vtx(I,1) = plaque1.vtx(I,1) - dl;

% Cut and translate plaque 2
I = ismember(plaque2.vtx,cut.vtx,'rows');
plaque2.vtx(I,1) = plaque2.vtx(I,1) + dl;

% Cut and translate Tjunction
I = ismember(Tjunction.vtx,int.vtx,'rows');
Tjunction.vtx(I,3) = Tjunction.vtx(I,3) + dl;

% Interfaces
I = ismember(int.vtx,cut.vtx,'rows');
int1 = int;
int1.vtx(I,1) = int1.vtx(I,1) - dl;
int2 = int;
int2.vtx(I,1) = int2.vtx(I,1) + dl;
int3 = int;
int3.vtx(:,3) = int3.vtx(:,3) + dl;

% Graphical representation
figure
% plot(mesh,'w')
hold on
% plot(bound,'r')
plot(plaque1,'g')
plot(plaque2,'r')
plot(Tjunction,'b')
plot(int1,'y')
plot(int2,'y')
plot(int3,'y')
axis equal
view(-50,20)

% Plaque with a holl
mesh = union(union(plaque1,plaque2),Tjunction);

% Dirichet condition
bound = mesh.bnd;
int   = union(union(int1,int2),int3);
dir   = setdiff(bound,int);

% Graphical representation
figure
plot(mesh,'w')
hold on
plot(int1,'g')
plot(int2,'r')
plot(int3,'b')
plot(dir,'y')
axis equal
view(-50,20)

% Domain
omega = dom(mesh,7);

% Finites elements space
u = fem(mesh,'P2');

% Dirichlet condition
u = dirichlet(u,dir);

% Junction
u = junction(u,int1,1,int2,-1,int3,1);

% Mass matrix
tic
M = integral(omega,u,u);
toc

% Rigidity matrix
tic
K = integral(omega,grad(u),grad(u));
toc

% Find eigen values
tic
[V,EV] = eigs(K,M,2*Neig,'SM');
toc

% Normalization
V = V./(max(max(abs(V))));

% Sort by ascending order
[EV,ind] = sort(sqrt(real(diag(EV))));
V        = V(:,ind);

% Graphical representation
figure
for n = 1:9
    subplot(3,3,n)
    surf(u,V(:,n))
    title(['k = ',num2str(EV(n))])
    axis equal off
    colorbar
end

% Analytical solutions of eigenvalues for an arbitrary cube
ref = zeros(Neig^2,1);
l = 1;
for i = 1:Neig
    for j = 1:Neig
        ref(l) = pi*sqrt( (i/L(1))^2 + (j/L(2))^2 );
        l = l+1;
    end
end
ref = sort(ref);
ref = ref(1:Neig);

% Error
sol = EV(1:Neig);
[ref sol abs(sol-ref)./ref]


disp('~~> Michto gypsilab !')


