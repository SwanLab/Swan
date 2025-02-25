%% Calculate the coordinates necessary to plot the colored data
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [xdata,xdataz,ydata,ydataz] = prepare_coordinates(elementType,xmesh,ymesh,doZ)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  Calculate positions where to plot and where the interpolation data is given.
%
%  P1 (DoF=3) - Element Numbering
%
%            3
%           /|
%          / |
%         /  |
%        /   |
%       /    |
%      /     |
%     /      |
%    /       |
%   1--------2
%
%  For P1 elements xdata=xdataz=xmesh and ydata=ydataz=ymesh.
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
%  1.) {1,2,4},{2,3,4},{1,4,3}
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
%      /     |
%     /      |
%    /       |
%   1---6----2
%
%  Depending on refinement the triangle {1,2,3} is divided into following
%  subtrinagles:
%
%  1.) {1,16,15},{16,62,c},{62,2,24},{24,43,c},{43,3,35},{35,15,c},
%      {15,1,16},{15,16,c},{62,24,c},{43,35,c} where c denotes the center
%      point of the main triangle {1,2,3}
%  2.) {1,6,5},{6,4,5},{6,2,4},{4,3,5}
%

%% Input Parameters
%
%  elementType: Type of the Lagrangian Finite Element (P1,P3,...)
%  xmesh:       holds the "nodal mesh" coordinates as given per
%               "savemesh" command.
%  ymesh:       holds the "nodal mesh" coordinates as given per
%               "savemesh" command.
%  doZ:         Create output for z-style plotting
%

%% Output Parameters
%
%  xdata:       the points given to the OpenGl driver to plot the PDE problem
%  xdataz:      the points to create a z-style meshplot
%  ydata:       the points given to the OpenGl driver to plot the PDE problem
%  ydataz;      the points to create a z-style meshplot
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
function [xdata,xdataz,ydata,ydataz] = prepare_coordinates(elementType,xmesh,ymesh,doZ)

    switch (elementType)
        case {'P0','P1'}
            xdata=xmesh;
            ydata=ymesh;
            %z-ploting for P1 is simply plotting at the (xmesh, ymesh) points
            xdataz=[];
            ydataz=[];
        case ('P1b')
            %Calculate the barycenter
            px4=(xmesh(3,:)+xmesh(2,:)+xmesh(1,:))/3;
            py4=(ymesh(3,:)+ymesh(2,:)+ymesh(1,:))/3;
            %Refine Mesh
            xmesh124=[xmesh(1,:);xmesh(2,:);px4]; ymesh124=[ymesh(1,:);ymesh(2,:);py4];
            xmesh234=[xmesh(2,:);xmesh(3,:);px4]; ymesh234=[ymesh(2,:);ymesh(3,:);py4];
            xmesh143=[xmesh(1,:);px4;xmesh(3,:)]; ymesh143=[ymesh(1,:);py4;ymesh(3,:)];
            %Assemble the refined mesh
            xdata=[xmesh124 xmesh234 xmesh143];
            ydata=[ymesh124 ymesh234 ymesh143];
            if doZ
                %z-ploting for P1b is plotting the triangle vertices = (xmesh, ymesh) points
                xdataz=xmesh;
                ydataz=ymesh;
            else
                xdataz=[];
                ydataz=[];
            end
        case ('P2')
            %Based on the nodal mesh points do refine the mesh according to the
            %extra points given in a P2 element
            if true
                %New point coordinates
                px15=xmesh(1,:)+(xmesh(3,:)-xmesh(1,:))/3;
                py15=ymesh(1,:)+(ymesh(3,:)-ymesh(1,:))/3;
                px35=xmesh(3,:)+(xmesh(1,:)-xmesh(3,:))/3;
                py35=ymesh(3,:)+(ymesh(1,:)-ymesh(3,:))/3;
                px43=xmesh(3,:)+(xmesh(2,:)-xmesh(3,:))/3;
                py43=ymesh(3,:)+(ymesh(2,:)-ymesh(3,:))/3;
                px24=xmesh(2,:)+(xmesh(3,:)-xmesh(2,:))/3;
                py24=ymesh(2,:)+(ymesh(3,:)-ymesh(2,:))/3;
                px62=xmesh(2,:)+(xmesh(1,:)-xmesh(2,:))/3;
                py62=ymesh(2,:)+(ymesh(1,:)-ymesh(2,:))/3;
                px16=xmesh(1,:)+(xmesh(2,:)-xmesh(1,:))/3;
                py16=ymesh(1,:)+(ymesh(2,:)-ymesh(1,:))/3;
                ctrx=(xmesh(1,:)+xmesh(2,:)+xmesh(3,:))/3;
                ctry=(ymesh(1,:)+ymesh(2,:)+ymesh(3,:))/3;
                %Refine the mesh
                xm1_16_15=[xmesh(1,:);px16;px15]; ym1_16_15=[ymesh(1,:);py16;py15];
                xm16_62_ctr=[px16;px62;ctrx]; ym16_62_ctr=[py16;py62;ctry];
                xm62_2_24=[px62;xmesh(2,:);px24]; ym62_2_24=[py62;ymesh(2,:);py24];
                xm24_43_ctr=[px24;px43;ctrx]; ym24_43_ctr=[py24;py43;ctry];
                xm43_3_35=[px43;xmesh(3,:);px35]; ym43_3_35=[py43;ymesh(3,:);py35];
                xm35_15_ctr=[px35;px15;ctrx]; ym35_15_ctr=[py35;py15;ctry];
                xm15_16_ctr=[px15;px16;ctrx]; ym15_16_ctr=[py15;py16;ctry];
                xm62_24_ctr=[px62;px24;ctrx]; ym62_24_ctr=[py62;py24;ctry];
                xm43_35_ctr=[px43;px35;ctrx]; ym43_35_ctr=[py43;py35;ctry];
                %Assemble the refined mesh
                xdata=[xm1_16_15,xm16_62_ctr,xm62_2_24,xm24_43_ctr,xm43_3_35, ...
                       xm35_15_ctr,xm15_16_ctr,xm62_24_ctr,xm43_35_ctr];
                ydata=[ym1_16_15,ym16_62_ctr,ym62_2_24,ym24_43_ctr,ym43_3_35, ...
                       ym35_15_ctr,ym15_16_ctr,ym62_24_ctr,ym43_35_ctr];
                if doZ
                    %For z-style plots we have to return a mesh of the triangle
                    %points={1,16,62,2,24,43,3,35,15}
                    xdataz=[xmesh(1,:);px16;px62;xmesh(2,:);px24;px43;xmesh(3,:);px35;px15];
                    ydataz=[ymesh(1,:);py16;py62;ymesh(2,:);py24;py43;ymesh(3,:);py35;py15];
                else
                    xdataz=[];
                    ydataz=[];
                end
            else
                px4=(xmesh(3,:)+xmesh(2,:))/2; py4=(ymesh(3,:)+ymesh(2,:))/2;
                px5=(xmesh(3,:)+xmesh(1,:))/2; py5=(ymesh(3,:)+ymesh(1,:))/2;
                px6=(xmesh(2,:)+xmesh(1,:))/2; py6=(ymesh(2,:)+ymesh(1,:))/2;
                xmesh165=[xmesh(1,:);px6;px5]; ymesh165=[ymesh(1,:);py6;py5];
                xmesh645=[px6;px4;px5]; ymesh645=[py6;py4;py5];
                xmesh624=[px6;xmesh(2,:);px4]; ymesh624=[py6;ymesh(2,:);py4];
                xmesh435=[px4;xmesh(3,:);px5]; ymesh435=[py4;ymesh(3,:);py5];
                xdata=[xmesh165 xmesh645 xmesh624 xmesh435];
                ydata=[ymesh165 ymesh645 ymesh624 ymesh435];
                if doZ
                    %Returns the Z-values at the triangle points {1,6,2,4,3,5}
                    xdataz=[xmesh(1,:);px6;xmesh(2,:);px4;xmesh(3,:);px5];
                    ydataz=[ymesh(1,:);py6;ymesh(2,:);py4;ymesh(3,:);py5];
                else
                    xdataz=[];
                    ydataz=[];
                end
            end
        otherwise
            error('Unknown Lagrangian Finite Element. Only P1 and P2 allowed');
    end

end