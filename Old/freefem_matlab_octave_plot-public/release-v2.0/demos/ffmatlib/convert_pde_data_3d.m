%% Detects the Lagrangian Finite Element Type and return the data in points tetrahedron format
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [elementType, xyzdata] = convert_pde_data_3d(points,tetrahedra,vhseq,xyzrawdata)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  Detect the type of the FE-Element and return the data in points
%  tetrahedron format (a horizontal matrix containing nT columns x 4 rows for P1
%  element data and nT columns x 10 rows for P2 element data).
%

%% Input Parameters
%
%  points:      A nodal coordinates points list created by savemesh()
%  tetrahedra:  A connectivity list created by savemesh()
%  xyzrawdata:  PDE FE-Space data exported by FreeFem++
%  vhseq:       Connectivity sequence of the FE-Space
%               (for non periodic P1 data this can be omitted)
%

%% Output Parameters
%
%  elementType: The type of the Lagrangian Finite Element {'P1','P2'}
%  xyzdata:     Data in points triangle format
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
function [elementType, xyzdata] = convert_pde_data_3d(points,tetrahedra,vhseq,xyzrawdata)

    [~,nt]=size(tetrahedra);
    [ndim,ndof]=size(xyzrawdata);
    [~,nv]=size(points);
    xyzdata=cell(1,ndim);
    %Detect Lagrangian Finite Element based on a check between vhseq and the number of triangles
    if (~isempty(vhseq))
        if (max(vhseq+1)~=length(xyzrawdata(1,:)))
            error('Datalength of PDE data ''XYData'' or ''FlowData'' does not match with Vh-Space connectivity ''VhSeq''');
        end
        switch (length(vhseq))
            case (nt)
                elementType='P0'; %nDoF=1 / element
                for i=1:ndim
                    cCols=xyzrawdata(i,:);
                    tmp=cCols(vhseq+1);
                    xyzdata{i}=[tmp;tmp;tmp;tmp];
                end
            case (4*nt)
                elementType='P1'; %nDoF=4 / element
                for i=1:ndim
                    cCols=xyzrawdata(i,:);
                    xyzdata{i}=reshape(cCols(vhseq+1),4,nt);
                end
            case (10*nt)
                elementType='P2'; %nDoF=10 / element
                for i=1:ndim
                    cCols=xyzrawdata(i,:);
                    xyzdata{i}=reshape(cCols(vhseq+1),10,nt);
                end
            otherwise
                error('Unknown Lagrangian Finite Element: ''VhSeq'' does not match with number of mesh triangles');
        end
    else
        %For P1 Data which comes not from a periodic problem there is a short cut without
        %vhseq because the order is already coded into the meshfile
        %Periodic problems must be processed with the vhseq option.
        if (ndof==nv)
            elementType='P1'; %nDoF=4 / element
            for i=1:ndim
                cCols=xyzrawdata(i,:);
                %short cut: P1 can data can be also converted without haven vhseq
                %vhseq already coded into the FE-Mesh
                xyzdata{i}=[cCols(tetrahedra(1,:)); cCols(tetrahedra(2,:)); cCols(tetrahedra(3,:)); cCols(tetrahedra(4,:))];
            end
        else
            error('Unknown FE-Space order: No ''VhSeq'' is given and nDoF does not match number of mesh vertices');
        end
    end

end
