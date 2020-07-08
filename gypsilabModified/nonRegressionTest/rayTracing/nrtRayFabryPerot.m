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
%|    #    |   FILE       : nrtRayFabryPerot.m                            |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray tracing with two planar reflecting meshes |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nvtx  = 1e2
Nray  = 1e5
Xsrc  = [0 0 0];
Xmes  = [0 0 0];
rad   = 0.5

% Fabry-Perot mesh
mesh     = mshSquare(Nvtx,[10 10]);
mesh.vtx = [-2+mesh.vtx(:,3) mesh.vtx(:,1) mesh.vtx(:,2)];
mesh2          = mesh;
mesh2.vtx(:,1) = mesh.vtx(:,1) + 4;
mesh           = union(mesh,mesh2);

% Graphical representation
[X,Y,Z] = sphere(50);
X = Xmes(1) + rad*X; Y = Xmes(2) + rad*Y; Z = Xmes(3) + rad*Z; 

% Plot mesh
figure
plot(mesh)
hold on
surf(X,Y,Z,zeros(size(X)));
axis equal;
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.3)

% Initialize ray tracing
ray = ray(mesh,Xsrc,Nray);
% plot(ray)

% Maximum distances 
r1000 = rad/2 * sqrt(Nray/1000);
r100  = rad/2 * sqrt(Nray/100);
r10   = rad/2 * sqrt(Nray/10);
rMax  = rad/2 * sqrt(Nray/2)

% Ray tracing
tic
ray = ray.tracer(100,rMax);
toc
% plot(ray)

% Spherical measures
tic
[I,src] = measure(ray,Xmes,rad,rMax);
toc

% Graphical representation for images
src = cell2mat(src);
plot3(src(:,1),src(:,2),src(:,3),'.')

% Energy
r   = 4*(1:length(I)-1);
ref = pi*rad^2./(4*pi*r.^2);
sol = zeros(size(ref));
for i = 1:length(ref)
    sol(i) = 0.5 * length(I{i+1})/Nray;
end

% Error
norm(ref-sol)/norm(ref)
norm(ref-sol,'inf')/norm(ref,'inf')

% Energy error in dB
figure
plot(r,10*log10(ref),'r',r,10*log10(sol),'+b')
hold on
plot([r1000 r1000],[-60,0],'k--')
text(r1000,-55,' n = 1000')
plot([r100 r100],[-60,0],'k--')
text(r100,-57,' n = 100')
plot([r10 r10],[-60,0],'k--')
text(r10,-55,' n = 10')
plot([rMax rMax],[-60,0],'k--')
text(rMax,-55,' n = 1')
xlabel('Distance source mesure')
ylabel('Energie mesuree (dB)')
title('Energie mesuree selon la distance de mesure')
legend({'Analytique','Statistique'})
grid on

% Erreur
figure
plot(r,10*log10(abs(sol-ref)./abs(ref)),'+-')
hold on
plot([r1000 r1000],[-60,-0],'k--')
text(r1000,-55,' n = 1000')
plot([r100 r100],[-60,-0],'k--')
text(r100,-57,' n = 100')
plot([r10 r10],[-60,-0],'k--')
text(r10,-55,' n = 10')
plot([rMax rMax],[-60,-0],'k--')
text(rMax,-55,' n = 1')
xlabel('Distance source mesure')
ylabel('Energie mesuree (dB)')
title('Energie mesuree selon la distance de mesure')
legend({'Analytique','Statistique'})
grid on


disp('~~> Michto gypsilab !')



