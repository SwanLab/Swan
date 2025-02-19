%movie_matlab.m Creates a movie
%
% Note: This Code is not compatible with Octave!
%
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

rgb=[0,0,0; 0,0,255; 0,255,255; 0, 255, 0; 255, 255, 0; 255, 0, 0; 224,0,255]/255;
[sz1,~]=size(rgb);
usedmap=interp1(linspace(0,1,sz1),rgb,linspace(0,1,200));

[p,b,t,nv,nbe,nt,labels] = ffreadmesh('movie.msh');
[qh]=ffreaddata('movie_qh.txt');
n=250;
for j=1:n
    name = sprintf('movie_temp_%i.txt', j+10000);
    u(:,j) = ffreaddata(name);
    name = sprintf('movie_psi_%i.txt', j+10000);
    psi(:,j) = ffreaddata(name);
end
%capture frames from the simulation data
figure('Resize','off','ToolBar','none','MenuBar','none','Position',[50 50 1000 230]);
set(gcf,'Renderer','OpenGL');
opengl('software'); %a getframe issue on older Matlab versions
hold on;
ffpdeplot(p,b,t,'VhSeq',qh,'XYData',u(:,1),'ColorMap',usedmap,'ColorRange', [0 50], ...
          'Title','Frame 0','CBTitle','dT[K]');
ffpdeplot(p,b,t,'VhSeq',qh,'XYData',psi(:,1),'XYStyle','off', ...
          'Contour','on','CGridParam',[250, 250],'CStyle','dashedneg', ...
          'ColorMap','off','ColorRange','off','ColorBar','off');
xlabel('x');
ylabel('y');
axis tight manual;
hax = gca;
%set(hax, 'NextPlot', 'replaceChildren');
m = moviein(n);
for j = 1:n
    cla;
    titlestr = sprintf('Free Convection Cavity [%0.1fsec]', (j-1)*0.04);
    ffpdeplot(p,b,t,'VhSeq',qh,'XYData',u(:,j),'ColorMap',usedmap,'ColorRange', [0 50], ...
          'Title',titlestr,'CBTitle','dT[K]');
    ffpdeplot(p,b,t,'VhSeq',qh,'XYData',psi(:,j),'XYStyle','off', ...
          'Contour','on','CGridParam',[250, 250],'CStyle','dashedneg', ...
          'ColorMap','off','ColorRange','off','ColorBar','off');
    drawnow;
    m(:,j) = getframe(gcf);
end
%create AVI file from animation
v = VideoWriter('movie.avi');
set(v, 'FrameRate', 15);
set(v, 'Quality', 95);
open(v);
for j = 1:n
    writeVideo(v, m(:,j));
end
close(v);
close all;
