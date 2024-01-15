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

[p,b,t,nv,nbe,nt,labels,regions]=ffreadmesh('region.msh');
[vh]=ffreaddata('region_vh.txt');
[u]=ffreaddata('region_data.txt');

figure;
ffpdeplot(p,b,t,'Mesh','on','Boundary','on');

figure;
ffpdeplot(p,b,t,'Mesh','on','Boundary','off','RLabels',regions(1));

plotcols = rand(numel(regions),3);
figure;
ffpdeplot(p,b,t,'Mesh','on','Boundary','off','RLabels',regions,'RColors',plotcols);

figure;
ffpdeplot(p,b,t,'VhSeq',vh,'XYData',u,'ZStyle','continuous','Mesh','on','MColor','k');
