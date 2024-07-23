%ffcplxmesh.m Returns a complex cartesian mesh grid
%
% Author: Chloros2 <chloros2@gmx.de>
% Created: 2018-11-11
%
%   This file is a part of the ffmatlib which is hosted at
%   https://github.com/samplemaker/freefem_matlab_octave_plot
%
%   [Z1, Z2] = ffcplxmesh(r1, r2, grid, numxy);
%
%   Creates a complex cartesian mesh grid spanned by the two complex
%   numbers r1, r2 and containing numxy grid lines
%
%   Parameters
%   r1, r2:  ranges r1=[Re(from),Im(from)]; r2=[Re(to),Im(to)]
%   numxy:   number of grid lines [numx, numy]
%   ngrid:   multiplier for grid points [ngridx, ngridy]
%
%   Example:
%   [Z1, Z2] = ffcplxmesh([0,0], [2*pi(),2*pi()], [3,3], [5,5]);
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
function [Z1, Z2]=ffcplxmesh(r1, r2, ngrid, numxy)
    %Re = const; Im varies
    x1 = linspace(r1(1), r2(1), numxy(1));
    y1 = linspace(r1(2), r2(2), ngrid(1)*numxy(1));
    [X1, Y1] = meshgrid(x1, y1);
    Z1 = X1 + 1i*Y1;
    %Im = const; Re varies
    x2 = linspace(r1(2), r2(2), numxy(2));
    y2 = linspace(r2(1), r1(1), ngrid(2)*numxy(2));
    [Y2, X2] = meshgrid(x2, y2);
    Z2 = X2 + 1i*Y2;
end
