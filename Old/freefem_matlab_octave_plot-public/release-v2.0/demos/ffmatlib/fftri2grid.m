%% Interpolates from a 2D triangular mesh to a 2D cartesian or curved grid
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [w1, ...] = fftri2grid (x, y, tx, ty, tu1, ...)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  fftri2grid computes the function values w1, w2, ... over a mesh grid defined
%  by the arguments x, y from a set of functions u1, u2, ... with values
%  given on a triangular mesh tx, ty. The values are computed using first order
%  or second order approximating basis functions (P1, P1b or P2 - Lagrangian Finite
%  Elements). The function values w1, w2, ... are real if tu1, tu2, ... are real
%  or complex if tu1, tu2, ... are complex. The mesh grid x, y can be cartesian
%  or curved. fftri2grid returns NaNs when an interpolation point is outside the
%  triangular mesh. fftri2gridfast.c is a runtime optimized mex implementation
%  of the function fftri2grid.m.
%

%% Triangle numbering
%
%  P1 (DoF=3) - Element Numbering
%
%            3
%           /|
%          / |
%         /  |
%        /   |
%       /    |
%      /  p  |
%     /      |
%    /       |
%   1--------2
%
%  For P1 elements coloring data is located at the mesh points.
%
%
%  P1b (DoF=4) - Element Numbering
%
%            3
%           /|
%          / |
%         /  |
%        /   |
%       /    |
%      /     |
%     /   4  |
%    /       |
%   1--------2
%
%  Note: The barycenter 4 is defined by the bubble function.
%
%  The triangle {1,2,3} is divided into following subtrinagles:
%
%
%  P2 (DoF=6) - Element Numbering
%
%            3
%           /|
%          / |
%         /  |
%        /   |
%       5    4
%      /  p  |
%     /      |
%    /       |
%   1---6----2
%

%% P1-Element Approximation
%
% The values inside the triangle are choosen by the use of the
% Ansatzfunktion:
%
% $$ u_p(x,y)=a_1+a_2x+a_3y $$
%
% Choose a coordinate transformation $$ (x,y) \rightarrow (w_1,w_2) $$ where
%
% $$ w_1=\frac{A_1(x,y)}{A_0} $$; $$ w_2=\frac{A_2(x,y)}{A_0} $$
%
% $$ A_1 $$ is the subtriangle area {p,2,3} and $$ A_2 $$ is the subtriangle area {p,3,1}
%
% $$ \rightarrow u_p(w_1,w_2)=c_1+c_2w_1+c_3w_2 $$
%
% This equation can be solved for the constants depending on $$ u $$ given at
% the triangle vertices. The base functions turn out to be:
%
% $$ N_1=w_1 $$
%
% $$ N_2=w_2 $$
%
% $$ N_3=(1-w_1-w_2) $$
%
% The interpolation of $$ u_p $$ at the point p with the coordinates $$ (x,y) $$ is
% calculated by the values given on the triangle edge/vertices by:
%
% $$ u_p=u_1N_1+u_2N_2+u_3N_3 $$
%

%% P1b-Element Approximation (P1 with bubble)
%
% Each triangle has one bubble (a point) at the triangle barycenter which
% adds a fourth degree of freedom. Properties of the bubble function:
%
%  - Vanishes on triangle edges
%  - Becomes 1 at the barycenter
%
% The base functions $$ N_1, N_2, N_3 $$ were already discussed in the
% P1-Element approximation above. The bubble function now becomes:
%
% $$ \beta(w_1,w_2) = 3^3 N_1 N_2 N_3 $$
%
% To interpolate inside a triangle the P1-Element interpolation
% is augmented with the bubble function. Since FreeFem++ returns the
% actual PDE value at the triangle vertices ($$ u_1,u_2,u_3 $$) and the
% barycenter ($$ u_4 $$) the interpolation of $$ u_p $$ at the point
% p with the coordinates $$ (x,y) $$ becomes:
%
% $$ u_p=u_1N_1+u_2N_2+u_3N_3+(u_4-\frac{1}{3}(u_1+u_2+u_3))\beta $$
%
% which yields:
%
% $$ u_p=u_1N_1+u_2N_2+u_3N_3+9(3u_4-(u_1+u_2+u_3))N_1 N_2 N_3 $$
%

%% P2-Element Approximation
%
% The values inside the triangle are choosen by the use of the
% Ansatzfunktion:
%
% $$ u_p(x,y)=a_1+a_2x+a_3y+a_4x^2+a_5xy+a_6y^2 $$
%
% Choose a coordinate transformation $$ (x,y) \rightarrow (w_1,w_2) $$ where
%
% $$ w_1=\frac{A_1(x,y)}{A_0} $$; $$ w_2=\frac{A_2(x,y)}{A_0} $$
%
% $$ A_1 $$ is the subtriangle area {p,2,3} and $$ A_2 $$ is the subtriangle area {p,3,1}
%
% $$ \rightarrow u_p(w_1,w_2)=c_1+c_2w_1+c_3w_2+c_4w_1^2+c_5w_1w_2+c_6w_2^2 $$
%
% This equation can be solved for the constants depending on $$ u $$ given at
% the triangle edge/vertices. The base functions turn out to be:
%
% $$ N_1=w_1(2w_1-1) $$
%
% $$ N_2=w_2(2w_2-1) $$
%
% $$ N_3=(1-w_1-w_2)(1-2w_1-2w_2) $$
%
% $$ N_4=4w_2(1-w_1-w_2) $$
%
% $$ N_5=4w_1(1-w_1-w_2) $$
%
% $$ N_6=4w_1w_2 $$
%
% The interpolation of $$ u_p $$ at the point p with the coordinates $$ (x,y) $$ is
% calculated by the values given on the triangle edge/vertices by:
%
% $$ u_p=u_1N_1+u_2N_2+u_3N_3+u_4N_4+u_5N_5+u_6N_6 $$
%

%% Input Parameters
%
%  x,y:       Meshgrid where to interpolate
%  tx,ty:     Triangular Mesh
%  varargin:  FE-Space data (Lagrangian Elements Type P1 or P2)

%% Output Parameters
%
%  varargout: Interpolation of the FE-Space data on the grid-points x,y
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
function [varargout] = fftri2grid(x, y, tx, ty, varargin)

    if nargin < 5
        error('Wrong number arguments');
    end
    %Number of variables passed to the interpolation function
    nDim=length(varargin);
    if ~isequal(size(x), size(y))
        error('Meshgrid x,y sizes must be equal');
    end
    if ~isequal(size(tx), size(ty))
        error('Nodal coordinate arrays tx,ty must be equal');
    end
    [npts,~]=size(tx);
    if (npts ~= 3)
        error('Nodal coordinate arrays must have 3 rows');
    end
    [ndof,~]=size(varargin{1}); %Autodetect Type of Lagrangian Finite Element
    [ny,nx]=size(x);
    varargout=cell(1,nDim);
    for i=1:nDim
        varargout{i}=NaN(ny,nx);
    end
    x1=tx(1,:);
    y1=ty(1,:);
    x2=tx(2,:);
    y2=ty(2,:);
    x3=tx(3,:);
    y3=ty(3,:);
    invA0=(1.0)./((y2-y3).*(x1-x3)+(x3-x2).*(y1-y3));
    for mx=1:nx
        for my=1:ny
            px=x(my,mx);
            py=y(my,mx);
            A1=((y2-y3).*(px-x3)+(x3-x2).*(py-y3)).*invA0;
            A2=((y3-y1).*(px-x3)+(x1-x3).*(py-y3)).*invA0;
            A3=1.0-A1-A2;
            pos=find(((A1>=-1e-13)&(A2>=-1e-13)&(A3>=-1e-13)),1,'first');
            if ~isempty(pos)
                w1=A1(pos);
                w2=A2(pos);
                w3=A3(pos);
                switch (ndof)
                    case 1 %P0 - Peace wise constant
                        for i=1:nDim
                            varargout{i}(my,mx)=varargin{i}(1,pos);
                        end
                    case 3 %P1 - Lagrangian Elements
                        for i=1:nDim
                            varargout{i}(my,mx)=varargin{i}(1,pos).*w1+ ...
                                                varargin{i}(2,pos).*w2+ ...
                                                varargin{i}(3,pos).*w3;
                        end
                    case 4 %P1b - (P1 with bubble)
                        for i=1:nDim
                            varargout{i}(my,mx)=varargin{i}(1,pos).*w1+ ...
                                                varargin{i}(2,pos).*w2+ ...
                                                varargin{i}(3,pos).*w3+ ...
                                                9*(3*varargin{i}(4,pos)- ...
                                                (varargin{i}(1,pos)+varargin{i}(2,pos)+varargin{i}(3,pos))).* ...
                                                w1.*w2.*w3;
                        end
                    case 6 %P2 - Lagrangian Elements
                        for i=1:nDim
                            varargout{i}(my,mx)=varargin{i}(1,pos).*w1.*(2*w1-1)+ ...
                                                varargin{i}(2,pos).*w2.*(2*w2-1)+ ...
                                                varargin{i}(3,pos).*(1-w1-w2).*(1-2*(w1+w2))+ ...
                                                varargin{i}(4,pos).*4.*w2.*(1-w1-w2)+ ...
                                                varargin{i}(5,pos).*4.*w1.*(1-w1-w2)+ ...
                                                varargin{i}(6,pos).*4.*w1.*w2;
                        end
                    %otherwise
                        %poor man's fall through
                end
            end
        end
    end

end
