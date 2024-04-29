%% Convert FreeFem++ mesh connectivity to nodal coordinates
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [xmesh,ymesh,zmesh] = prepare_mesh_3d(points,tetrahedra)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Licence
%
%  Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
%  This program is free software: you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see
%  <https://www.gnu.org/licenses/>.
%

%% Code
function [xmesh,ymesh,zmesh] = prepare_mesh_3d(points,tetrahedra)
    xpts=points(1,:);
    ypts=points(2,:);
    zpts=points(3,:);
    xmesh=[xpts(tetrahedra(1,:)); xpts(tetrahedra(2,:)); xpts(tetrahedra(3,:)); xpts(tetrahedra(4,:))];
    ymesh=[ypts(tetrahedra(1,:)); ypts(tetrahedra(2,:)); ypts(tetrahedra(3,:)); ypts(tetrahedra(4,:))];
    zmesh=[zpts(tetrahedra(1,:)); zpts(tetrahedra(2,:)); zpts(tetrahedra(3,:)); zpts(tetrahedra(4,:))];
end
