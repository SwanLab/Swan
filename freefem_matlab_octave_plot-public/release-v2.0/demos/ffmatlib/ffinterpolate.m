%% Interpolates PDE simulation data on a cartesian or curved meshgrid
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [w1, ...] = ffinterpolate(p, ~, t, vhseq, x, y, u1, ...)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  ffinterpolate interpolates the real or complex valued multidimensional data
%  u1, ... given on a triangular mesh defined by the points p, triangles t and
%  boundary b, onto a cartesian or curved grid defined by the arguments x, y.
%  The return values w1, ... are real if u1, ... are real or complex
%  if u1, ... are complex. ffinterpolate() is a wrapper function that
%  calls the library functions fftri2grid() and fftri2gridfast().
%

%% Input Parameters
%
%  x,y:       Meshgrid where to interpolate
%  p:         A nodal coordinates points list created by savemesh()
%  t:         A connectivity list created by savemesh()
%  varargin:  FE-Space data (multidimensional vectorfield)
%  vhseq:     Sequence of the FE-Space Function in terms of connectivity


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
function [varargout] = ffinterpolate(p, ~, t, vhseq, x, y, varargin)

    if nargin < 7
        error('Wrong number arguments');
    end

    %\TODO: Better solution?
    nDim=length(varargin);
    [~,len]=size(rowvec(varargin{1}));
    xyrawdata=zeros(nDim,len);
    for i=1:nDim
        xyrawdata(i,:)=rowvec(varargin{i});
    end
    [~,pdeData]=convert_pde_data(p,t,vhseq,xyrawdata);
    [xmesh,~,ymesh,~]=prepare_mesh(p,t);

    if exist('fftri2gridfast','file')
        [varargout{1:nargout}] = fftri2gridfast(x,y,xmesh,ymesh,pdeData{:});
    else %[a,b]=udata{:}
        fprintf('Note: To improve runtime build MEX function fftri2gridfast() from fftri2gridfast.c\n');
        [varargout{1:nargout}] = fftri2grid(x,y,xmesh,ymesh,pdeData{:});
    end

 end

 function [S] = rowvec(S)
    [sz1,sz2]=size(S);
    if sz1>sz2
        S=S.';
    end
end
