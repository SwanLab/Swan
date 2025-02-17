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
%|    #    |   FILE       : nrtRayCube.m                                  |
%|    #    |   VERSION    : 0.41                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 01.04.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray tracing with absorbing cube               |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
L     = [5 4 3]
% Xsrc  = [3.9 2.8 1.7]
% Xmes  = [2.1 2.2 1.3]
Xsrc  = [4 2 1.7];
Xmes  = [2 2 1.7];
Nray  = 1e5
rad   = 0.2

% Read mesh
mesh = msh('cube1e1.mesh');

% Mesh reshape to feat dimensions in L
mesh.vtx = 0.5 .* (1 + mesh.vtx);
mesh.vtx = (ones(size(mesh.vtx,1),1)*L) .* mesh.vtx;

% Initialize ray
ray = ray(mesh,Xsrc,Nray);

% Graphical sphere
[X,Y,Z] = sphere(50);
X = Xmes(1) + rad*X; Y = Xmes(2) + rad*Y; Z = Xmes(3) + rad*Z; 

% Graphical representation
figure
plot(mesh)
hold on
% plot(ray)
surf(X,Y,Z,ones(size(X)))
axis equal;
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.3)
colorbar

% Material properties (http://www.odeon.dk/material-manufactures)
ray.msh.col(mesh.col==1) = 10;    % X = 1
ray.msh.col(mesh.col==2) = 20;    % Y = 1
ray.msh.col(mesh.col==3) = 30;    % X = 0
ray.msh.col(mesh.col==4) = 40;    % Y = 0
ray.msh.col(mesh.col==5) = 50;    % Z = 0
ray.msh.col(mesh.col==6) = 60;    % Z = 1
% ray.msh.col(:) = 60;

% Maximum distances 
r1000 = rad/2 * sqrt(Nray/1000);
r100  = rad/2 * sqrt(Nray/100);
r10   = rad/2 * sqrt(Nray/10);
rMax  = rad/2 * sqrt(Nray/2);

% Ray-tracing
tic
ray = ray.tracer(100,rMax);
toc
% plot(ray)

% Images sources
tic
[img,nrg] = ray.image(Xmes,rad,rMax);
toc

% Graphical representation for images
n   = size(img,1);
tmp = ones(n,1)*Xmes + img;
tmp = msh(tmp,(1:n)');
plot(tmp,10*log10(nrg(:,1)));
axis equal
colorbar

% Analytical solution
load('odeon.mat')
num = [30 ; % X = 0
    10 ;    % X = 1
    40 ;    % Y = 0
    20 ;    % Y = 1
    50 ;    % Z = 0
    60];    % Z = 1
% num = 60 * ones(6,1);
mat = odeon.mat(num,:);
rhm = 30;
air = 5.5 * (50/rhm) .* (odeon.frq/1000).^1.7 * 1e-4;
tic
[imgRef,nrgRef] = rayCubeAnalytic(L,mat,air,Xsrc,Xmes,rMax);
toc

% Graphical representation
plot3(imgRef(:,1)+Xmes(1),imgRef(:,2)+Xmes(2),imgRef(:,3)+Xmes(3),'ok')

% Data
r    = sqrt(sum(img.^2,2));
rRef = sqrt(sum(imgRef.^2,2));
sol  = mean(nrg,2);
ref  = mean(nrgRef,2);

% Energy error in dB
figure
plot(rRef,10*log10(ref),'or',r,10*log10(sol),'+b')
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

% Audio file
[audio,fs] = audioread('anechoicVoice.wav');

% Fir from images
T = floor(r/340*fs) + 1;
for i = 1:size(nrg,2)
    fir8(:,i) = accumarray(T,nrg(:,i),[max(T) 1]);
end

% Bank filtering
dir8 = ray.bank(256,fs);
rir  = 0;
abs(sum(sum(dir8,2)) - 1)
for i = 1:size(dir8,2)
    rir = rir + fftfilt(dir8(:,i),fir8(:,i));
%     figure
%     freqz(dir8(:,i),1,2*256)
end

% Graphical representation
figure
plot(rir)

% Audio rendering (5 seconds)
out = fftfilt(rir,audio(1:5*fs));
% sound(out,fs)



disp('~~> Michto gypsilab !')


