% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique  
% Algiane Froehly, CARDAMOME, INRIA-SOFT 
% LGPL Lesser General Public License v3.0. 
% Remeshing using Mmg tools : https://www.mmgtools.org      

% Cleaning
clear all
close all
clc

% Gypsilab path
%run('../../addpathGypsilab.m')

% Parameters
Nvtx = 1e4;            % Number initial of vertices
x0   = [0 0];          % Initial source position
L    = [3 3];          % Room size
tf   = 2;              % Final time
c    = 1;              % Sound celerity
emin = 5e-2;           % Minimum edge size

% Initial condition for magnitude and speed
ui = @(X) exp(-(X(:,1)-x0(1)).^2/5e-2) .* exp(-(X(:,2)-x0(2)).^2/5e-2);  
vi = @(X) zeros(size(X,1),1);                                            

% Time step (CFL)
dx = 0.1%emin;          
dt = 0.25*dx(1)/c;
dt = tf/(length(0:dt:tf)+1);

% Distance function
dst = @(u) emin./(abs(u)/max(abs(u))+emin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ORIGINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial geometry and fem objects
mesh  = mshSquare(Nvtx,L);
u     = fem(mesh,'P1');

% Initial magnitude and speed (u0 and v0)
u0 = ui(u.unk);
v0 = vi(u.unk);
u1 = u0 + dt * v0;

% Graphical representation
figure(1)
subplot(1,2,1)
plot(mesh,'w'); 
hold on
plot(mesh,u1);
alpha(0.5)
colorbar
% caxis([0 0.5])
axis equal;
title('Solution')
xlabel('X');   ylabel('Y');   zlabel('Z');

subplot(1,2,2)
plot(mesh,'w'); 
hold on
plot(mesh,dst(u1));
alpha(0.5)
colorbar
% caxis([0 0.5])
axis equal;
title('Distance')
xlabel('X');   ylabel('Y');   zlabel('Z');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REFINMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mmg refinment
refine = mmg(mesh);
hmin(refine,emin);
map(refine,dst(u1));
hausd(refine,1e-3)
mesh = run(refine);

% Geometry and fem objects
omega = dom(mesh,3);
u     = fem(mesh,'P1');

% Operators
M = integral(omega,u,u);
K = integral(omega,grad(u),grad(u));

% Initial magnitude and speed (u0 and v0)
u0 = ui(u.unk);
v0 = vi(u.unk);
u1 = u0 + dt * v0;

% Graphical representation
figure(2)
subplot(1,2,1)
plot(mesh,'w'); 
hold on
plot(mesh,u1);
alpha(0.5)
colorbar
caxis([0 0.5])
axis equal;
title('Solution')
xlabel('X');   ylabel('Y');   zlabel('Z');

subplot(1,2,2)
plot(mesh,'w'); 
hold on
plot(mesh,dst(u1));
alpha(0.5)
colorbar
% caxis([0 0.5])
axis equal;
title('Distance')
xlabel('X');   ylabel('Y');   zlabel('Z');

% % Video output
% vidObj = VideoWriter('1Ddirichlet.avi');
% open(vidObj);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order two sheme in times
t = dt;
while t < tf
    tic
    
    % Leap-frog
    u2 = (M+(c*dt)^2*K) \ (M*(2*u1-u0));
    t  = t + dt;
    
    toc
    tic
    
    % Graphical representation
    figure(2)
    clf
    subplot(1,2,1)
    plot(mesh,'w');
    hold on
    plot(mesh,u2);
    alpha(0.5)
    colorbar
    caxis([0 0.2])
    axis equal;
    title('Solution')
    xlabel('X');   ylabel('Y');   zlabel('Z');
    
    subplot(1,2,2)
    plot(mesh,'w');
    hold on
    plot(mesh,dst(u2));
    alpha(0.5)
    colorbar
    % caxis([0 0.5])
    axis equal;
    title('Distance')
    xlabel('X');   ylabel('Y');   zlabel('Z');

    %     currFrame = getframe(gcf);
    %     writeVideo(vidObj,currFrame);    
    
    toc
    tic

    % Mmg refinment
    refine = mmg(mesh);
    hmin(refine,emin);
    map(refine,dst(u2));
    verbose(refine,-1)
    new = run(refine);
    
    toc
    tic
    
    % Data interpolation
    X  = u.unk;
    F0 = scatteredInterpolant(X(:,1),X(:,2),u1);
    u0 = F0(new.vtx(:,1),new.vtx(:,2));%   feval(u,u1,mesh);
    F1 = scatteredInterpolant(X(:,1),X(:,2),u2);
    u1 = F1(new.vtx(:,1),new.vtx(:,2));%   feval(u,u1,mesh);
    
    toc
    tic
    
    % Change mesh and fem
    mesh  = new;
    omega = dom(mesh,3);
    u     = fem(mesh,'P1');
    
    % Change operators
    M = integral(omega,u,u);
    K = integral(omega,grad(u),grad(u));
    
    toc
    '================================================================='
end

%  % Close the file.
%  close(vidObj);


disp('~~> Michto gypsilab !')



