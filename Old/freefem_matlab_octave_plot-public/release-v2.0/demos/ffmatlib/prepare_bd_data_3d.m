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
%Convert the FE-Spacefunction data at the boundaries to euclidean coordinates
%Data in points/vertex format (P1-Simulation): Works for P1 FE-space only
function [elementType, varargout] = prepare_bd_data_3d(vhseq,nv,triangles,tetrahedra,xyzrawdata,bdlabels)
    [~,nt]=size(tetrahedra);
    [ndim,ndof]=size(xyzrawdata);
    if (~isempty(vhseq))
        if (max(vhseq+1)~=length(xyzrawdata(1,:)))
            error('Datalength of PDE data ''XYData'' or ''FlowData'' does not match with Vh-Space connectivity ''VhSeq''');
        end
        switch (length(vhseq))
            case (nt)
                elementType='P0'; %nDoF=1 / element
            case (4*nt)
                elementType='P1'; %nDoF=4 / element
            case (10*nt)
                elementType='P2'; %nDoF=10 / element
            otherwise
                error('Unknown Lagrangian Finite Element: ''VhSeq'' does not match with number of mesh triangles');
        end
    else
        if (ndof==nv)
            elementType='P1'; %nDoF=4 / element
        else
            error('Unknown FE-Space order: No ''VhSeq'' is given and nDoF does not match number of mesh vertices');
        end
    end

    varargout=cell(1,ndim);

    switch (elementType)
        case ('P0')
            %Hack: Create lookup table
            %NOTE: Displays incorrectly at the borders
            [~,nt]=size(tetrahedra);
            ptnum=tetrahedra(1:4,:);
            tmp=[vhseq+1,vhseq+1,vhseq+1,vhseq+1]';
            lookup=zeros(nv,1);
            lookup(ptnum(:))=tmp(:);
            %Convert only the boundaries specified in bdlabels
            if ~isempty(bdlabels)
                keep=(triangles(4,:)==bdlabels(1));
                for i=2:numel(bdlabels)
                    keep=(keep | (triangles(4,:)==bdlabels(i)));
                end
                for i=1:ndim
                    cols=xyzrawdata(i,:);
                    varargout{i}=[cols(lookup(triangles(1,keep))); cols(lookup(triangles(2,keep))); cols(lookup(triangles(3,keep)))];
                end
            else
                for i=1:ndim
                    cols=xyzrawdata(i,:);
                    varargout{i}=[cols(lookup(triangles(1,:))); cols(lookup(triangles(2,:))); cols(lookup(triangles(3,:)))];
                end
            end
        case ('P1')
            if (ndof==nv)
                %Convert only the boundaries specified in bdlabels
                if ~isempty(bdlabels)
                    keep=(triangles(4,:)==bdlabels(1));
                    for i=2:numel(bdlabels)
                        keep=(keep | (triangles(4,:)==bdlabels(i)));
                    end
                    for i=1:ndim
                        cols=xyzrawdata(i,:);
                        varargout{i}=[cols(triangles(1,keep)); cols(triangles(2,keep)); cols(triangles(3,keep))];
                    end
                else
                    for i=1:ndim
                        cols=xyzrawdata(i,:);
                        varargout{i}=[cols(triangles(1,:)); cols(triangles(2,:)); cols(triangles(3,:))];
                    end
                end
            else
                error('boundary: unable to recognize input data format');
            end
        case ('P2')
            %Hack: Create lookup table
            %Note: Display surface only as P1 (extract P1 subdata from P2)
            [~,nt]=size(tetrahedra); %mesh numbering + region in terms of connectivity
            ptnum=tetrahedra(1:4,:); %remove region information, keep vertice mesh data
            tmp=reshape(vhseq+1,10,nt); %bring into triangle format
            vheven=tmp(1:4,:); %get the vertice corner points, remove higher order rubbish
            lookup=zeros(nv,1);
            lookup(ptnum(:))=vheven(:); %create index where to find a specific vertice number
            %Convert only the boundaries specified in bdlabels
            if ~isempty(bdlabels)
                keep=(triangles(4,:)==bdlabels(1));
                for i=2:numel(bdlabels)
                    keep=(keep | (triangles(4,:)==bdlabels(i)));
                end
                for i=1:ndim
                    cols=xyzrawdata(i,:);
                    varargout{i}=[cols(lookup(triangles(1,keep))); cols(lookup(triangles(2,keep))); cols(lookup(triangles(3,keep)))];
                end
            else
                for i=1:ndim
                    cols=xyzrawdata(i,:);
                    varargout{i}=[cols(lookup(triangles(1,:))); cols(lookup(triangles(2,:))); cols(lookup(triangles(3,:)))];
                end
            end
        otherwise
            error('Unknown Lagrangian Finite Element. Only P0, P1 and P2 allowed');
    end
end