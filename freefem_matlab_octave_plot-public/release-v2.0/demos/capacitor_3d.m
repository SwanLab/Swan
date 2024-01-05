% 3D parallel plate capacitor problem
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

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('cap3d.mesh');
[vh]=ffreaddata('cap3vh.txt');
[u]=ffreaddata('cap3dpot.txt');
[Ex,Ey,Ez]=ffreaddata('cap3dvec.txt');

%Definition in edp file
CK=30;
CA=31;
 
S1=[-0 0.375 0.0; ...
    0.375 0 0.0];
S2=[0.0 0.375 0.5; ...
    0.375 0 0.5];
S3=[0.75 0.375 0.0; ...
    0.375 0.75 0.0];

figure;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,'Slice',S1,S2,S3,'SGridParam',[50,50], ...
            'Boundary','on','BDLabels',[CK,CA],'XYZStyle','monochrome', ...
            'ColorMap',jet(200),'ColorBar','on','BoundingBox','on');
ylabel('y');
xlabel('x');
zlabel('z');

S1=[-0 0.0 0.0];
S2=[0.0 0.375 0.5];
S3=[0.75 0.375 0.0];

figure;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,'Slice',S1,S2,S3,'SGridParam',[50,50], 'Project2D', 'on', ...
            'Boundary','off','ColorMap',jet(200),'ColorBar','on');

figure;
ffpdeplot3D(p,b,t,'VhSeq',vh,'FlowData',[Ex,Ey,Ez],'FGridParam3D',[8,8,5], ...
            'FLim3D',[0.125,0.625;0.125,0.625;0.1,0.4],'BDLabels',[CK,CA], ...
            'XYZStyle','monochrome');

figure;
ffpdeplot3D(p,b,t,'XYZStyle','monochrome');
str={sprintf('nVertex: %i',nv);
     sprintf('nTets:   %i',nt);
     sprintf('nTris:   %i',nbe)};
annotation('textbox',[0.05 0.05 0.2 0.15],'String',str, ...
           'FitBoxToText','on','FontName','Courier','EdgeColor',[1 1 1]);
