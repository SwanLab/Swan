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
%|    #    |   FILE       : nrtFemWave2D.m                                |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2018                                    |
%| ( === ) |   SYNOPSIS   : 2D wave propagation                           |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nvtx = 1e3;            % Number of vertices
x0   = [0 0];          % Initial position
L    = [1 1];          % Strinf size
tf   = 2;              % Final time
c    = 1;              % Sound celerity

% Initial condition for magnitude and speed
ui = @(X) exp(-(X(:,1)-x0(1)).^2/1e-1) .* exp(-(X(:,2)-x0(2)).^2/1e-1);  
vi = @(X) zeros(size(X,1),1);                                            

% Meshes
% mesh = mshSquare(Nvtx,L);
mesh = mshDisk(Nvtx,1);

% Boundary
bound = mesh.bnd;

% Domain
omega = dom(mesh,3);

% Finites elements space
u = fem(mesh,'P1');

% Dirichlet
% u = dirichlet(u,bound);

% Graphical representation
plot(mesh,'w'); 
hold on
plot(bound,'r')
plot(omega)
plot(u,'go')
hold off
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.99)

% Time step (CFL)
dx = mesh.stp;          
dt = 0.25*dx(1)/c;
dt = tf/(length(0:dt:tf)+1);

% Mass matrix
M = integral(omega,u,u);
sum(sum(M)) - pi

% Rigidity matrix
K = integral(omega,grad(u),grad(u));
sum(sum(K))

% Initial magnitude and speed (u0 and v0)
u0 = ui(u.unk);
v0 = vi(u.unk);

% First step scheme (Forward Euler) 
u1 = u0 + dt * v0;

% Graphical representation
figure(10)
graph(u,u1);
colorbar
axis equal
caxis([0 0.5])
view(0,60)
camlight
material dull

% % Video output
% vidObj = VideoWriter('1Ddirichlet.avi');
% open(vidObj);

% Order two sheme in times
t = dt;
while t < tf
    % Leap-frog
    u2 = M \ ( M*(2*u1-u0) - (c*dt)^2*K*u1 );
    t  = t + dt;
    
    % Graphical representation
    figure(10)
    clf
    graph(u,u2);
    colorbar
    axis equal
    caxis([0 0.5])
    view(0,60)
    camlight
    material dull
    title(['Wave propagation t = ',num2str(t,'%6.3f')])
    %     currFrame = getframe(gcf);
    %     writeVideo(vidObj,currFrame);
    pause(dt)
    
    % Incrementation
    u0 = u1;
    u1 = u2;  
end

%  % Close the file.
%  close(vidObj);


disp('~~> Michto gypsilab !')



