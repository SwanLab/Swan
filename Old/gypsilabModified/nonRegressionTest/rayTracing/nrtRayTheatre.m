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
%|    #    |   FILE       : nrtRayTheatre.m                               |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray tracing with theatre                      |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
Xsrc  = [0 8 1.7];
Xmes  = [0 -20 8];
Nray  = 1e5
rad   = 1

% Read mesh
mesh        = msh('theatre.mesh');
mesh.col(:) = 100; % Rough concrete (Bobran, 1973)

% Extract toit
ctr = mesh.ctr;
mesh = mesh.sub(ctr(:,3)<14);

% Initialize ray
ray1 = ray(mesh,Xsrc,Nray);

% Graphical sphere
[X,Y,Z] = sphere(50);
X = Xmes(1) + rad*X; Y = Xmes(2) + rad*Y; Z = Xmes(3) + rad*Z; 

% Graphical representation
figure
plot(mesh,'w')
hold on
plot(ray1)
surf(X,Y,Z,ones(size(X)))
axis equal
xlabel('X');   ylabel('Y');   zlabel('Z');
view(0,45)
alpha(0.5)

% Maximum distances 
r1000 = rad/2 * sqrt(Nray/1000);
r100  = rad/2 * sqrt(Nray/100);
r10   = rad/2 * sqrt(Nray/10);
rMax  = rad/2 * sqrt(Nray/2);

% Ray-tracing
tic
ray1 = ray1.tracer(30,rMax);
toc
% plot(ray1)

% Images sources
tic
[img,nrg] = ray1.image(Xmes,rad,rMax);
toc

% Data
r   = sqrt(sum(img.^2,2));
sol = mean(nrg,2);

% Impacts
ray2 = ray(mesh,Xmes,img);
ray2 = ray2.tracer;
plot(ray2)
plot(msh(ray2.pos{2},(1:length(ray2))',10*log10(sol)));
axis equal
colorbar

% Energy in dB
figure
plot(r,10*log10(sol),'+b')
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
grid on

% Audio file
[audio,fs] = audioread('anechoicVoice.wav');

% Fir from images
T = floor(r/340*fs) + 1;
for i = 1:size(nrg,2)
    fir8(:,i) = accumarray(T,nrg(:,i),[max(T) 1]);
end

% Bank filtering
dir8 = ray1.bank(256,fs);
rir  = 0;
for i = 1:size(dir8,2)
    rir = rir + fftfilt(dir8(:,i),fir8(:,i));
end

% Graphical representation
figure
plot(rir)

% Audio rendering (5 seconds)
out = fftfilt(rir,audio(1:5*fs));
% sound(out,fs)



disp('~~> Michto gypsilab !')


