% periodic boundary condition problem
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-10-31
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

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('periodic.msh');
[vh]=ffreaddata('periodic_vh.txt');
[u]=ffreaddata('periodic.txt');

%%%%%% Mesh Plot

figure;

ffpdeplot(p,b,t,'Mesh','on','Boundary','on','Title','Mesh');

ylabel('y');
xlabel('x');

%%%%%% Show Labels

figure;

ffpdeplot(p,b,t,'Mesh','off','Boundary','on','BDLabels',labels,'BDShowText','on', ...
                'Title','Boundary Labels');

%%%%%% 2D Patch (density map) Plot + Contour

figure;

ffpdeplot(p,b,t,'XYData',u,'VhSeq',vh,'Mesh','off', ...
          'Boundary','on','CBTitle','u','Contour','on', ...
          'Title','Patch Plot + Contour');

ylabel('y');
xlabel('x');