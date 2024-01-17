%% Interpolates from a 3D tetrahedral mesh to a 2D grid
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [w1, ...] = fftet2grid (x, y, z, tx, ty, tz, tu1, ...)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  fftet2grid computes the function values w1, w2, ... over a mesh grid defined
%  by the arguments x, y, z from a set of functions u1, u2, ... with values
%  given on a tetrahedral mesh tx, ty, tz. The values are computed using first order
%  or second order approximating basis functions (P1 or P2 - Lagrangian Finite
%  Elements). The function values w1, w2, ... are real if tu1, tu2, ... are real
%  or complex if tu1, tu2, ... are complex. The mesh grid x, y, z  can be cartesian
%  or curved. fftet2grid returns NaNs when an interpolation point is outside the
%  tetrahedral mesh. fftet2gridfast.c is a runtime optimized mex implementation
%  of the function fftet2grid.m.
%

%% Tetrahedron numbering
%
%  P1 (DoF=4) - Element Numbering
%
%               2
%              /\\
%             /  \ \
%            /    \  \
%           /      \   \
%          /        \    \
%         /          \     \
%        /            \      \
%       /      p       \       \ 3
%      /                \      /
%     /                  \    /
%    /                    \  /
%   /______________________\/
%  4                        1
%
%  For P1 elements coloring data is located at the mesh points.
%
%  P2 (DoF=10) - Element Numbering
%
%               2
%              /\\
%             /  \ \
%            /    \  \
%           /      \   \
%          /        \    8
%         /          \     \
%        9            5      \
%       /      p       \      _\ 3
%      /                \      /
%     /          -10-    \    6
%    /    _               \  /
%   /__________ ___________\/
%  4           7            1
%

%% P1-Element Approximation
%
% Choose a coordinate transformation $$ (x_1,x_2,x_3) \rightarrow (w_1,w_2,w_3) $$ where
%
% $$ w_1=\frac{V_1(x,y)}{V_0} $$; $$ w_2=\frac{V_2(x,y)}{V_0} $$; $$ w_3=\frac{V_3(x,y)}{V_0} $$
%
% $$ V_1 $$ is the subvolume {p,3,2,4}, $$ V_2 $$ the subvolume {p,4,1,3}, $$ V_3 $$ the subvolume {p,1,2,4}
%
% The values inside the tetrahedron are choosen by the use of the
% Ansatzfunktion:
%
% $$ \rightarrow u_p(w_1,w_2,w_3)=c_1+c_2w_1+c_3w_2+c_4w_3 $$
%
% This equation can be solved for the constants depending on $$ u $$ given at
% the triangle edge/vertices. The base functions turn out to be:
%
% $$ N_1=w_1(2w_1-1) $$
%
% $$ N_2=w_2(2w_2-1) $$
%
% $$ N_3=w_3(2w_3-1) $$
%
% $$ N_4=(1-w_1-w_2-w_3) $$
%
% The interpolation of $$ u_p $$ at the point p with the coordinates $$ (x,y,z) $$ is
% calculated by the values given on the tetrahedron edge/vertices by:
%
% $$ u_p=u_1N_1+u_2N_2+u_3N_3+u_4N_4 $$
%

%% P2-Element Approximation
%
% Choose a coordinate transformation $$ (x_1,x_2,x_3) \rightarrow (w_1,w_2,w_3) $$ where
%
% $$ w_1=\frac{V_1(x,y)}{V_0} $$; $$ w_2=\frac{V_2(x,y)}{V_0} $$; $$ w_3=\frac{V_3(x,y)}{V_0} $$
%
% $$ V_1 $$ is the subvolume {p,3,2,4}, $$ V_2 $$ the subvolume {p,4,1,3}, $$ V_3 $$ the subvolume {p,1,2,4}
%
% The values inside the tetrahedron are choosen by the use of the
% Ansatzfunktion:
%
% $$ \rightarrow u_p(w_1,w_2,w_3)=c_1+c_2w_1+c_3w_2+c_4w_3+c_5w_1w_2+c_6w_1w_3+c_7w_2w_3+c_8w_1^2+c_9w_2^2+c_{10}w_3^2 $$
%
% This equation can be solved for the constants depending on $$ u $$ given at
% the tetrahedron edge/vertices. The base functions turn out to be:
%
% $$ N_1=w_1(2w_1-1) $$
%
% $$ N_2=w_2(2w_2-1) $$
%
% $$ N_3=w_3(2w_3-1) $$
%
% $$ N_4=(1-w_1-w_2-w_3)(1-2w_1-2w_2-2w_3) $$
%
% $$ N_5=4w_1w_2 $$
%
% $$ N_6=4w_1w_3 $$
%
% $$N_7=4w_1(1-w_1-w_2-w_3)$$
%
% $$ N_8=4w_2w_3 $$
%
% $$ N_9=4w_2(1-w_1-w_2-w_3) $$
%
% $$ N_{10}=4w_3(1-w_1-w_2-w_3) $$
%
% The interpolation of $$ u_p $$ at the point p with the coordinates $$ (x,y,z) $$ is
% calculated by the values given on the tetrahedron edge/vertices by:
%
% $$ u_p=u_1N_1+u_2N_2+u_3N_3+u_4N_4+u_5N_5+u_6N_6+u_7N_7+u_8N_8+u_9N_9+u_{10}N_{10} $$
%

%% Input Parameters
%
%  x,y:       Meshgrid where to interpolate
%  tx,ty:     Triangular Mesh
%  varargin:  FE-Space data (Lagrangian Elements P1 or P2)

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
function [varargout] = fftet2grid(x, y, z, tx, ty, tz, varargin)

    if nargin < 7
        error('Wrong number arguments');
    end
    %Number of variables passed to the interpolation function
    nDim=length(varargin);
    if ~isequal(size(x), size(y), size(z))
        error('Meshgrid x,y,z sizes must be equal');
    end
    if ~isequal(size(tx), size(ty), size(tz))
        error('Nodal coordinate arrays tx,ty,tz must be equal');
    end
    [npts,ntet]=size(tx);
    if (npts ~= 4)
        error('Nodal coordinate arrays must have 4 rows');
    end
    [ndof,~]=size(varargin{1}); %Autodetect Type of Lagrangian Finite Element
    [ny,nx]=size(x);
    varargout=cell(1,nDim);
    for i=1:nDim
        varargout{i}=NaN(ny,nx);
    end
    p1=[tx(1,:);ty(1,:);tz(1,:)];
    p2=[tx(2,:);ty(2,:);tz(2,:)];
    p3=[tx(3,:);ty(3,:);tz(3,:)];
    p4=[tx(4,:);ty(4,:);tz(4,:)];
    invVTET=abs(1.0./(dot(cross(p3-p1,p2-p1),p4-p1,1)));
    %Convert tetrahedrons into enclosing cubes
    txmx=max(max(tx(1,:),tx(2,:)),max(tx(3,:),tx(4,:)));
    txmn=min(min(tx(1,:),tx(2,:)),min(tx(3,:),tx(4,:)));
    tymx=max(max(ty(1,:),ty(2,:)),max(ty(3,:),ty(4,:)));
    tymn=min(min(ty(1,:),ty(2,:)),min(ty(3,:),ty(4,:)));
    tzmx=max(max(tz(1,:),tz(2,:)),max(tz(3,:),tz(4,:)));
    tzmn=min(min(tz(1,:),tz(2,:)),min(tz(3,:),tz(4,:)));
    for mx=1:nx
        for my=1:ny
            preselect=((x(my,mx) <= txmx) & ...
                       (x(my,mx) >= txmn) & ...
                       (y(my,mx) <= tymx) & ...
                       (y(my,mx) >= tymn) & ...
                       (z(my,mx) <= tzmx) & ...
                       (z(my,mx) >= tzmn));
            pp1=p1(:,preselect); pp2=p2(:,preselect);
            pp3=p3(:,preselect); pp4=p4(:,preselect);
            invVTETT=invVTET(:,preselect);
            [~,sz2]=size(invVTETT);
            O=ones(1,sz2);
            Xp=[x(my,mx)*O; y(my,mx)*O; z(my,mx)*O];
            V1=volume(pp4-pp2,pp3-pp2,Xp-pp2).*invVTETT;
            V2=volume(pp3-pp1,pp4-pp1,Xp-pp1).*invVTETT;
            V3=volume(pp4-pp1,pp2-pp1,Xp-pp1).*invVTETT;
            %V1=dot(cross(pp4-pp2,pp3-pp2),Xp-pp2,1).*invVTETT;
            %V2=dot(cross(pp3-pp1,pp4-pp1),Xp-pp1,1).*invVTETT;
            %V3=dot(cross(pp4-pp1,pp2-pp1),Xp-pp1,1).*invVTETT;
            V4=1.0-V1-V2-V3;
            pos=find(((V1>=-1e-13) & (V2>=-1e-13) & (V3>=-1e-13) & (V4>=-1e-13)),1,'first');
            if ~isempty(pos)
              for i=1:nDim
                U{i}=varargin{i}(:,preselect);
              end
                w1=V1(pos);
                w2=V2(pos);
                w3=V3(pos);
                w4=V4(pos);
                switch (ndof)
                    case 1 %P0 - Peacewise constant
                        for i=1:nDim
                            varargout{i}(my,mx)=(U{i}(1,pos)+ ...
                                                 U{i}(2,pos)+ ...
                                                 U{i}(3,pos)+ ...
                                                 U{i}(4,pos))/4.0;
                        end
                    case 4 %P1 - Lagrangian Elements
                        for i=1:nDim
                            varargout{i}(my,mx)=U{i}(1,pos).*w1+ ...
                                                U{i}(2,pos).*w2+ ...
                                                U{i}(3,pos).*w3+ ...
                                                U{i}(4,pos).*w4;
                        end
                    case 10 %P2 - Lagrangian Elements
                        for i=1:nDim
                            varargout{i}(my,mx)=U{i}(1,pos).*w1.*(2*w1-1)+ ...
                                                U{i}(2,pos).*w2.*(2*w2-1)+ ...
                                                U{i}(3,pos).*w3.*(2*w3-1)+ ...
                                                U{i}(4,pos).*(1-w1-w2-w3).*(1-2*(w1+w2+w3))+ ...
                                                U{i}(5,pos).*4*w1.*w2+ ...
                                                U{i}(6,pos).*4*w1.*w3+ ...
                                                U{i}(7,pos).*4*w1.*(1-w1-w2-w3)+ ...
                                                U{i}(8,pos).*4*w2.*w3+ ...
                                                U{i}(9,pos).*4*w2.*(1-w1-w2-w3)+ ...
                                                U{i}(10,pos).*4*w3.*(1-w1-w2-w3);
                        end
                    %otherwise
                        %poor man's fall through
                end
            end
        end
    end

end

function [vol]=volume(A,B,C)
   AxB=[A(2,:).*B(3,:)-A(3,:).*B(2,:);A(3,:).*B(1,:)-B(3,:).*A(1,:);A(1,:).*B(2,:)-B(1,:).*A(2,:)];
   D=AxB.*C;
   vol=D(1,:)+D(2,:)+D(3,:);
end
