% Create a movie or anim gif of the heat transfer in a von Karman Vortex Street
% Runtime: ~25min
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2019-03-25
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

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('karman_vortex.msh');

[w]=ffreaddata('karman_vortex_vorticity_10000.txt');
[u]=ffreaddata('karman_vortex_temperature_10000.txt');
[pp]=ffreaddata('karman_vortex_pressure_10000.txt');
[xh]=ffreaddata('karman_vortex_xh.txt');
[mh]=ffreaddata('karman_vortex_mh.txt');

rgb=[0,0,0; 0,0,255; 0,255,255; 0, 255, 0; 255, 255, 0; 255, 0, 0; 224,0,255]/255;
[sz1,~]=size(rgb);
usedmap1=interp1(linspace(0,1,sz1),rgb,linspace(0,1,200));

rgb=[255,0,0; 255,255,255; 0,0,255]/255;
[sz1,~]=size(rgb);
usedmap2=interp1(linspace(0,1,sz1),rgb,linspace(0,1,200));

fig=figure('Position', [20 50 1000 700]);
annotation(fig,'textbox',[0.46 0.04 0.7 0.02], ...
           'String','https://github.com/samplemaker/freefem_matlab_octave_plot/', ...
           'FitBoxToText','off','Color','k','FontName','Courier','FontWeight','normal', ...
           'FontSize',9,'EdgeColor','none','Interpreter', 'none');

hax1=subplot(3,1,1);

ffpdeplot(p,b,t,'VhSeq',xh,'XYData',w, 'Mesh','off', ...
          'ColorRange',[-0.1,0.1],'ColorMap',usedmap1, ...
          'Boundary','off','CBTitle','w','Title','Vorticity');
axis tight;

hax2=subplot(3,1,2);
ffpdeplot(p,b,t,'VhSeq',xh,'XYData',u,'Mesh','off', ...
          'ColorRange',[0.0,0.3],'ColorMap',usedmap1, ...
          'Boundary','off','CBTitle','dT','Title','Temperature');
axis tight;

hax3=subplot(3,1,3);
ffpdeplot(p,b,t,'VhSeq',mh,'XYData',pp,'Mesh','off', ...
          'ColorRange', [-0.005,0.005],'ColorMap',usedmap2, ...
          'Boundary','off','CBTitle','P','Title','Pressure');
axis tight;

n=699;
for j = 1:n
    w_name = sprintf('karman_vortex_vorticity_%i.txt', j+10000-1);
    u_name = sprintf('karman_vortex_temperature_%i.txt', j+10000-1);
    p_name = sprintf('karman_vortex_pressure_%i.txt', j+10000-1);
    [w]=ffreaddata(w_name);
    [u]=ffreaddata(u_name);
    [pp]=ffreaddata(p_name);
    set(fig,'CurrentAxes',hax1);
    ffpdeplot(p,b,t,'VhSeq',xh,'XYData',w, 'Mesh','off', ...
              'ColorRange',[-0.1,0.1],'ColorMap',usedmap1, ...
              'Boundary','off','CBTitle','w','Title','Vorticity');
    set(fig,'CurrentAxes',hax2);
    ffpdeplot(p,b,t,'VhSeq',xh,'XYData',u,'Mesh','off', ...
              'ColorRange',[0.0,0.3],'ColorMap',usedmap1, ...
              'Boundary','off','CBTitle','dT','Title','Temperature');
    set(fig,'CurrentAxes',hax3);
    ffpdeplot(p,b,t,'VhSeq',mh,'XYData',pp,'Mesh','off', ...
              'ColorRange', [-0.005,0.005],'ColorMap',usedmap2, ...
              'Boundary','off','CBTitle','P','Title','Pressure');

    drawnow;

    %Export each frame to png
    file = sprintf('movie_%d.png', 10000+j-1);
    print(gcf,'-dpng','-r100',file);
    M(j) = getframe(fig);
end

%At this point you can create an avi file but then you have
%to run this file in Matlab
%Or you can use Octave and Image Magic to create an animated
%give from the *.png files issueing:
%convert -delay 5x150 *.png anim.gif

v = VideoWriter('movie.avi');
set(v, 'FrameRate', 15);
set(v, 'Quality', 95);
open(v);
for j = 1:n
    writeVideo(v, M(j));
end
close(v);
close all;

