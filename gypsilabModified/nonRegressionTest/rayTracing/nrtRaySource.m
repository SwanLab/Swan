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
%|    #    |   FILE       : nrtRaySource.m                                |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray tracing with spherical measure            |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Nray = 1e5
Xsrc = [0 0 0]
Xmes = [1 0 0]
rad  = 0.5

% Initialize ray tracing
ray = ray([],Xsrc,Nray);

% Spherical measurement
tic
I = inSphere(ray,Xmes,rad);
toc
Ip    = true(Nray,1);
Ip(I) = 0;

% Graphical representation
[X,Y,Z] = sphere(50);
X = Xmes(1) + rad*X; Y = Xmes(2) + rad*Y; Z = Xmes(3) + rad*Z; 

% Plot raytracing
figure
plot(ray.sub(Ip))
hold on
surf(X,Y,Z,zeros(size(X)));
quiver3(ray.pos{end}(I,1),ray.pos{end}(I,2),ray.pos{end}(I,3),...
    ray.dir(I,1),ray.dir(I,2),ray.dir(I,3),'g')
axis equal;
xlabel('X');   ylabel('Y');   zlabel('Z');
grid on
alpha(0.3)

% Maximum distances 
r1000 = rad/2 * sqrt(Nray/1000);
r100  = rad/2 * sqrt(Nray/100);
r10   = rad/2 * sqrt(Nray/10);
rMax  = rad/2 * sqrt(Nray/2);

% Energy computation for a range of distances
r = 4*rad:(rMax-4*rad)/100:rMax;
n = zeros(size(r));
E = zeros(size(r));
for i = 1:length(r)
    % Distance
    tmp  = [r(i) 0 0];
    
    % Measures
    I = ray.inSphere(tmp,rad);
    
    % Number of ray measured
    n(i) = length(I);
    
    % Statistical energy
    E(i) = n(i)/Nray;    
end

% Analytical energy
ref = pi*rad^2./(4*pi*r.^2);

% Energy error in dB
figure
plot(r,10*log10(ref),'r',r,10*log10(E),'+b')
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
plot(r,10*log10(abs(E-ref)./abs(ref)),'+-')
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
ylabel('Erreur relative (dB)')
title('Energie mesuree selon la distance de mesure')
grid on


disp('~~> Michto gypsilab !')

