% A heat transfer problem
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2019-03-01
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

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('heat_transfer.msh');
[vh]=ffreaddata('heat_transfer_vh.txt');
[u]=ffreaddata('heat_transfer_data.txt');

%%%%%% 3D Surf Plot Gouraud lighting

figure;

ffpdeplot(p,b,t, ...
          'VhSeq',vh, ...
          'XYData',u, ...
          'ZStyle','continuous', ...
          'ColorMap',jet(150), ...
          'Mesh','off', ...
          'CBTitle','U', ...
          'Title','Gouraud');
ylabel('y');
xlabel('x');
zlabel('u');
grid;

lighting gouraud;
view([-47,24]);
camlight('headlight');

axis tight;
