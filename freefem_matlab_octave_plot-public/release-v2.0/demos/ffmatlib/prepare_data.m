%% Convert raw FE-Space data to plot- and interpolation data
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [varargout] = prepare_data (elementType,triangles,xydata,xmesh,ymesh,z)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  prepare_data calculates the PDE data matrices which are used to create the
%  patch color plots. The color data can be a scalar or a vector field,
%  depending on the input (number of columns) xydata.
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
%  triangles:   Connectivity list of the triangular mesh
%  xydata:      FE-Space data (can be multidimensional [ux,uy])
%  xmesh,ymesh: Nodal coordinates of the mesh
%  doZ:         Create output for z-style plotting
%

%% Output Parameters
%
%  varargout: Data/interpolation of the FE-Space data on
%  1.) the points given to the OpenGl driver to plot the PDE problem
%  2.) the points to create a z-style meshplot
%  3.) the points used for further interpolation code (for P2 a matrix
%  containing 6 rows x nTriangle columns and for P1 3rows x nTriangle
%  columns
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
function [varargout] = prepare_data(elementType,triangles,xydata,xmesh,ymesh,doZ)

    ndim=numel(xydata);
    %Three outputs:
    %1. Refined mesh data
    %2. Refined mesh for z-plot data
    %3. Data to be used for fftri2grid() to create contour plots
    varargout=cell(1,3*ndim);
    switch (elementType)
        case {'P0','P1'}
            %For P1 elements interpolation data and zdata is equal to the
            %mesh data.
            varargout=xydata;
            for i=1:ndim
                %z-style plots
                varargout{i+ndim}=[];
                %interpolation code / fftri2grid()
                varargout{i+2*ndim}=[];
            end
        case ('P1b')
            for i=1:ndim
               %cc is a nTriangle columns x 4 rows matrix containing
               %all FE-Space data (= 3 triangle vertices + 1x bubble
               %function)
               cc=xydata{i};
               %the value at the barycenter is the bubblevalue
               uintrp=cc(4,:);
               %Refine the mesh
               um1_2_4=[cc(1,:);cc(2,:);uintrp];
               um2_3_4=[cc(2,:);cc(3,:);uintrp];
               um1_4_3=[cc(1,:);uintrp;cc(3,:)];
               %Assemble the refined mesh
               cdata=[um1_2_4,um2_3_4,um1_4_3];
               varargout{i}=cdata;
               %z-ploting for P1b is plotting the triangle vertices = (xmesh, ymesh) points
               varargout{i+ndim}=cc(1:3,:);
               %used by fftri2grid()
               varargout{i+2*ndim}=cc;
            end
        %TODO: Remove "else" part and improve speed by moving all xmesh-calculations before the for loop!!!
        case ('P2')
            [~,nt]=size(triangles);
            for i=1:ndim
                if true
                    %cc is a nTriangle columns x 6 rows matrix containing
                    %all FE-Space data
                    cc=xydata{i};
                    ax=xmesh(1,:);bx=xmesh(2,:);cx=xmesh(3,:);
                    ay=ymesh(1,:);by=ymesh(2,:);cy=ymesh(3,:);
                    invA0=(1.0)./((by-cy).*(ax-cx)+(cx-bx).*(ay-cy));
                    px(1,:)=xmesh(1,:)+(xmesh(3,:)-xmesh(1,:))/3;
                    py(1,:)=ymesh(1,:)+(ymesh(3,:)-ymesh(1,:))/3; %p15
                    px(2,:)=xmesh(3,:)+(xmesh(1,:)-xmesh(3,:))/3;
                    py(2,:)=ymesh(3,:)+(ymesh(1,:)-ymesh(3,:))/3; %p35
                    px(3,:)=xmesh(3,:)+(xmesh(2,:)-xmesh(3,:))/3;
                    py(3,:)=ymesh(3,:)+(ymesh(2,:)-ymesh(3,:))/3; %p43
                    px(4,:)=xmesh(2,:)+(xmesh(3,:)-xmesh(2,:))/3;
                    py(4,:)=ymesh(2,:)+(ymesh(3,:)-ymesh(2,:))/3; %p24
                    px(5,:)=xmesh(2,:)+(xmesh(1,:)-xmesh(2,:))/3;
                    py(5,:)=ymesh(2,:)+(ymesh(1,:)-ymesh(2,:))/3; %p62
                    px(6,:)=xmesh(1,:)+(xmesh(2,:)-xmesh(1,:))/3;
                    py(6,:)=ymesh(1,:)+(ymesh(2,:)-ymesh(1,:))/3; %p16
                    px(7,:)=(xmesh(1,:)+xmesh(2,:)+xmesh(3,:))/3;
                    py(7,:)=(ymesh(1,:)+ymesh(2,:)+ymesh(3,:))/3; %center
                    uintrp=zeros(7,nt);
                    for j=1:7
                        w1=((by-cy).*(px(j,:)-cx)+(cx-bx).*(py(j,:)-cy)).*invA0;
                        w2=((cy-ay).*(px(j,:)-cx)+(ax-cx).*(py(j,:)-cy)).*invA0;
                        uintrp(j,:)=cc(1,:).*w1.*(2*w1-1)+cc(2,:).*w2.*(2*w2-1)+ ...
                                    cc(3,:).*(1-w1-w2).*(1-2*(w1+w2))+ ...
                                    cc(4,:).*4.*w2.*(1-w1-w2)+ ...
                                    cc(5,:).*4.*w1.*(1-w1-w2)+cc(6,:).*4.*w1.*w2;
                    end
                    %Refine the mesh
                    um1_16_15=[cc(1,:);uintrp(6,:);uintrp(1,:)];
                    um16_62_ctr=[uintrp(6,:);uintrp(5,:);uintrp(7,:)];
                    um62_2_24=[uintrp(5,:);cc(2,:);uintrp(4,:)];
                    um24_43_ctr=[uintrp(4,:);uintrp(3,:);uintrp(7,:)];
                    um43_3_35=[uintrp(3,:);cc(3,:);uintrp(2,:)];
                    um35_15_ctr=[uintrp(2,:);uintrp(1,:);uintrp(7,:)];
                    um15_16_ctr=[uintrp(1,:);uintrp(6,:);uintrp(7,:)];
                    um62_24_ctr=[uintrp(5,:);uintrp(4,:);uintrp(7,:)];
                    um43_35_ctr=[uintrp(3,:);uintrp(2,:);uintrp(7,:)];
                    %Assemble the refined mesh
                    cdata=[um1_16_15,um16_62_ctr,um62_2_24,um24_43_ctr,um43_3_35, ...
                           um35_15_ctr,um15_16_ctr,um62_24_ctr,um43_35_ctr];
                    varargout{i}=cdata;
                    %z-style plots. Returns the z-values (=color values) at the triangle
                    %points {1,6,2,4,3,5}
                    if doZ
                        varargout{i+ndim}=[cc(1,:);uintrp(6,:);uintrp(5,:);cc(2,:); ...
                                           uintrp(4,:);uintrp(3,:);cc(3,:);uintrp(2,:);uintrp(1,:)];
                    else
                        varargout{i+ndim}=[];
                    end
                    %used by fftri2grid()
                    varargout{i+2*ndim}=cc;
                else
                    %c=reshape(xydata{i},3,nt*2);
                    %The corner points of the main triangles {1,2,3}
                    %cEven=c(:,1:2:nt*2);
                    %The corner points of the subtriangles {4,5,6}
                    %cOdd=c(:,2:2:nt*2);
                    cEven=xydata{i}(1:3,:);
                    cOdd=xydata{i}(4:6,:);
                    c165=[cEven(1,:);cOdd((6-3),:);cOdd((5-3),:)];
                    c645=[cOdd((6-3),:);cOdd((4-3),:);cOdd((5-3),:)];
                    c624=[cOdd((6-3),:);cEven(2,:);cOdd((4-3),:)];
                    c435=[cOdd((4-3),:);cEven(3,:);cOdd((5-3),:)];
                    varargout{i}=[c165 c645 c624 c435];
                    if doZ
                        %Returns the z-values (=color values) at the triangle
                        %points {1,6,2,4,3,5}
                        varargout{i+ndim}=[cEven(1,:);cOdd((6-3),:);cEven(2,:); ...
                                           cOdd((4-3),:);cEven(3,:);cOdd((5-3),:)];
                    else
                        varargout{i+ndim}=[];
                    end
                    %used by fftri2grid()
                    varargout{i+2*ndim}=xydata{i};
                end
            end
    otherwise
        error('Unknown Lagrangian Finite Element. Only P1, P1b and P2 allowed');
    end

end

