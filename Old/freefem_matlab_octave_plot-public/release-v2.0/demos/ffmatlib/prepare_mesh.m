%% Convert FreeFem++ mesh connectivity to nodal coordinates
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [xmesh,xpts,ymesh,ypts] = prepare_mesh(points,triangles)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  prepare_mesh creates triangle matrices from from a FreeFem++ mesh
%  created by the savemesh() command.
%

%% Input Parameters
%
%  points:    A nodal coordinates points list created by savemesh()
%  triangles: A connectivity list created by savemesh()
%

%% Output Parameters
%
%  xmesh: 3xnTriangle matrix containing the x coordinates of the triangles
%  xpts:  The nodal coordinates
%  ymesh: 3xnTriangle matrix containing the x coordinates of the triangles
%  ypts:  The nodal coordinates
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
function [xmesh,xpts,ymesh,ypts] = prepare_mesh(points,triangles)
    xpts=points(1,:);
    ypts=points(2,:);
    xmesh=[xpts(triangles(1,:)); xpts(triangles(2,:)); xpts(triangles(3,:))];
    ymesh=[ypts(triangles(1,:)); ypts(triangles(2,:)); ypts(triangles(3,:))];
end
