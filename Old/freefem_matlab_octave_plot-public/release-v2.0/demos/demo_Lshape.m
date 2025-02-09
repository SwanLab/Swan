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
[p,b,t]=ffreadmesh('demo_Lshape_mesh.msh');
[vhseq]=ffreaddata('demo_Lshape_vh.txt');
[u,vx,vy]=ffreaddata('demo_Lshape_data.txt');

%%%%%%% 2D Mesh Plot
figure();
ffpdeplot(p,b,t, ...
          'Mesh','on', ...
          'Title','Mesh without boundary');

%%%%%%% Including Boundary
figure();
ffpdeplot(p,b,t, ...
          'Mesh','on', ...
          'Boundary','on', ...
          'Title','Mesh with boundary');

%%%%%%% Boundary only
figure();
ffpdeplot(p,b,t, ...
          'Boundary','on', ...
          'Title','Boundary only');

%%%%%%% 2D PDE Map Plot
figure();
ffpdeplot(p,b,t, ...
          'XYData',u, ...
          'VhSeq',vhseq, ...
          'Title','2D Density Plot');

%%%%%%% 2D PDE ColorMap Plot
figure();
mymap='jet';
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ColorMap',mymap, ...
          'Title',strcat('2D Density Plot - ColMap: ',mymap));

%%%%%%% 2D PDE Map Plot - customized HSV Map
figure();
hsv=[4./6.,1,0.5; 4./6.,1,1; 5./6.,1,1; 1,1,1; 1,0.5,1];
[sz1,~]=size(hsv);
cm_data=interp1(linspace(0,1,sz1),hsv,linspace(0,1,64));
usedmap=hsv2rgb(cm_data);
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq, ...
          'ColorMap',usedmap, ...
          'Title','2D Density Plot HSV custom MAP');

%%%%%%% 2D PDE Map Plot with Mesh
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Mesh','on', ...
          'Title','2D Density Plot with Mesh');

%%%%%%% 2D PDE Map Plot with Mesh
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Mesh','on', 'XYStyle', 'off', ...
          'ColorBar','off', 'ColorMap','off', 'ColorRange','off', ...
          'Title','2D Density Plot with Mesh without color (XYStyle)');

%%%%%%% 2D PDE Map Plot with Mesh without Colorbar
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Mesh','on','ColorBar','off', ...
          'Title','2D Density Plot with Mesh without Colbar');

%%%%%%% 2D PDE Map Plot with Mesh with Colorbar North
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Mesh','on','ColorBar','on', ...
          'ColorBar', 'northoutside', ...
          'Title','2D Density Plot with Mesh Colbar north');

%3D Surf Plots%

%%%%%%% 3D Surface
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ZStyle','continuous', ...
          'Title','Classic 3D Surf Plot');

%%%%%%% 3D Surface Lighting
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ZStyle','continuous', ...
          'Title','Classic 3D Surf Plot Gouraud');

lighting gouraud;
view([-166,40]);
camlight('left');
grid;

%%%%%%% 3D Surface Map-Patch including Mesh
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ZStyle','continuous','Mesh','on', ...
          'Title','3D Surf Plot including Mesh');

%%%%%%% 3D Mesh monochrome
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ZStyle','continuous','Mesh','on', 'XYStyle','off', ...
          'ColorBar','off', 'ColorMap','off', 'ColorRange','off', ...
          'Title','3D Mesh Plot - Monochrome (XYStyle)');

%%%%%%% 3D colored Mesh
% Not implemented

%Contour Plots%

%%%%%%% Contour Plot wo patch monochrome
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Boundary','on','Contour','on', 'XYStyle', 'off', 'CColor', [0,0,1], ...
          'ColorBar','off', 'ColorMap','off', 'ColorRange','off', ...
          'Title','Monochrome (blue) Contour without patch');

%%%%%%% Contour Plot wo patch flat (colored isolevels)
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Boundary','on','Contour','on', 'XYStyle', 'off', 'CColor', 'flat', ...
          'ColorMap','jet', ...
          'Title','Flat Contour without patch (colored isolevels)');

%%%%%%% Contour Plot Mixed with 2D Patch
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Boundary','on','Contour','on', ...
          'Title','Contour mixed with Patch Plot');

%%%%%%% Contour Plot Mixed with 2D Patch dashed lines
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Boundary','on','Contour','on', 'CStyle', 'dashed', ...
          'Title','Contour mixed with Patch Plot - dashed lines');

%%%%%%% Contour Plot
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'Boundary','on', ...
          'Contour','on','CGridParam',[10,10],'CLevels',15, ...
          'Title','Contour Plot with ugly set GridParam');

%%%%%%% Contour Plot with Labels
figure();
[handles,clab]=ffpdeplot(p,b,t, ...
                         'XYData',u,'VhSeq',vhseq,'Boundary','on','Contour','on','CLevels',5, ...
                         'Title','Contour mixed with Patch Plot with Labels');

texth=clabel(clab,handles(3),'fontsize', 8);
for i=1:size(texth)
    textstr=get(texth(i),'String');
    textnum=str2double(textstr);
    textstrnew=sprintf('%0.3f', textnum);
    set(texth(i),'String',textstrnew);
end

%%%%%%% 2D PDE Map Plot axis limits
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'XLim',[-0.25 0.75],'YLim',[0.25 1.25], ...
          'Title','2D Map Patch Plot with Limits');

%%%%%%% 3D Surface axis limits
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ZStyle','continuous', ...
          'XLim',[-0.25 0.75],'YLim',[0.25 1.25],'ZLim',[-0.1 0.1], ...
          'Title','3D Plot with xy and z-Limits');

%%%%%%% 3D Surface aspect ratio
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ZStyle','continuous','DAspect',[1 1 0.02], ...
          'Title','3D Plot Aspect ratio');

%%%%%%% 2D PDE Map Plot ColorRange
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ColorRange',[0 0.1], ...
          'Title','3D Plot Color Range set Test');

figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ColorRange','cropminmax', ...
          'XLim',[0 0.2],'YLim',[0.75 1.25], ...
          'Title','Autorange Cropped Section');

figure;
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ColorRange','cropcentered', ...
          'XLim',[0 0.2],'YLim',[0.75 1.25], ...
          'Title','Autorange Cropped Section, symmetric 0');

figure;
ffpdeplot(p,b,t, ...
          'XYData',u,'VhSeq',vhseq,'ColorRange','centered', ...
          'Title','Symmetric o');

%%%%%%% Quiver
figure();
ffpdeplot(p,b,t, ...
          'FlowData',[vx,vy],'VhSeq',vhseq,'Boundary','on', ...
          'Title','Quiver default');

axis tight;

%%%%%%% Quiver
figure();
ffpdeplot(p,b,t, ...
          'FlowData',[vx,vy],'VhSeq',vhseq,'FGridParam',[26,30],'Boundary','on', ...
          'Title','Quiver Plot with GridParam');

axis tight;

%%%%%%% Superposition Quiver + 2D
figure();
ffpdeplot(p,b,t, ...
          'XYData',u,'FlowData',[vx,vy],'VhSeq',vhseq,'FGridParam',[26,30],'Boundary','on', ...
          'Title','Quiver Plot combined with 2D Map Patch Plot');

axis tight;
