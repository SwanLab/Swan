% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-05-13
%
% Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.
%

clear all;
addpath('ffmatlib');


%Compute a conformal plot based on a PDE simulation result
[p,b,t]=ffreadmesh('complex_test.msh');
[u]=ffreaddata('complex_data.txt');
[vh]=ffreaddata('complex_vh.txt');
figure;
hold on;
%Z1: real part = const
%Z2: imag part = const
[Z1, Z2] = ffcplxmesh([0,0], [2*pi(),2*pi()], [3,3], [12,12]);
[w] = ffinterpolate(p,b,t,vh,real(Z1),imag(Z1),u);
plot(real(w),imag(w),'b','LineWidth',1);
[w] = ffinterpolate(p,b,t,vh,real(Z2),imag(Z2),u);
plot(real(w),imag(w),'r','LineWidth',1);
grid;
title('Conformal plot of a PDE solution');
daspect([1,1,1]);


%Project a polar plot onto a PDE solution
[p,b,t]=ffreadmesh('capacitor_2d.msh');
[u,~,~]=ffreaddata('capacitor_data_2d.txt');
[vh]=ffreaddata('capacitor_vh_2d.txt');
[Z1, Z2] = ffcplxmesh([0.5,0], [4,2*pi()], [10,10], [10,11]);
%Map to polar coordinates
strFunc='@(Z)(real(Z).*exp(1i*imag(Z)))';
f = str2func(strFunc);
U1 = f(Z1);
U2 = f(Z2);
figure;
hold on;
%constant radius
[w] = ffinterpolate(p,b,t,vh,real(U1),imag(U1),u);
plot3(real(U1),imag(U1),w,'g','LineWidth',1);
%constant angle
[w] = ffinterpolate(p,b,t,vh,real(U2),imag(U2),u);
plot3(real(U2),imag(U2),w,'g','LineWidth',1);
ffpdeplot(p,b,t, ...
          'VhSeq',vh, ...
          'XYData',u, ...
          'ZStyle','continuous', ...
          'Mesh','off', ...
          'CBTitle','U[V]', ...
          'Title','Curved Interpolation');
ylabel('y');
xlabel('x');
zlabel('u');
grid;
lighting gouraud;
view([-47,24]);
camlight('headlight');


%Interpolates on a curved line
N = 100;
s = linspace(0,2*pi(),N);
Z = 3.5*(cos(s)+1i*sin(s)).*sin(0.5*s);
w = ffinterpolate(p,b,t,vh,real(Z),imag(Z),u);
figure('position', [0, 0, 800, 300])
subplot(1,2,1);
hold on;
plot3(real(Z),imag(Z),real(w),'g','LineWidth',1);
ffpdeplot(p,b,t, ...
          'VhSeq',vh, ...
          'XYData',u, ...
          'ZStyle','continuous', ...
          'Mesh','off', ...
          'ColorBar','off', ...
          'Title','Single curve Interpolation');
ylabel('y');
xlabel('x');
zlabel('u');
grid;
lighting gouraud;
view([-47,24]);
camlight('headlight');
subplot(1,2,2);
plot(s,real(w),'b');
grid;
xlabel('s');
ylabel('u');
title('Interpolation Values');
