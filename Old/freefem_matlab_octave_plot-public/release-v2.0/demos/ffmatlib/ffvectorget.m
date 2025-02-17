%%  Extract from vectorized FE-Space
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2019-02-26
%
% [vhOut,uOut] = ffvectorget(struct,vhIn,u,elementNo)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  Gets component data and component FE-Space from a vectorized FE-Space
%

%% Input Parameters
%
%  elementNo: Index of component to be extracted
%  struct:    Structure of the vectorized FE-Space.
%             I.e. for "fespace Vh(Th, [P2,P1])" enter {'P2','P1'}
%  vhIn:      Vector valued input FE-Space
%  u:         Vector valued PDE solution
%

%% Output Parameters
%
%  vhOut: The component of the subspace/suborder
%  uOut:  The component of the PDE solution
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

%Extract the FE-Space for a particular element from a vectorized FE-Space
function [vhOut,uOut] = ffvectorget(struct,vhIn,u,elementNo)

    if (elementNo > numel(struct)) || (elementNo <= 0)
        error('Invalid element Index');
    end
 
    %DoFs of implemented vector spaces
    EnDoFs = {{'P0','P1','P1b','P2'},{1,3,4,6}};
    pattern = [];
    nDoF = [];
    %Cycle through all components and gather the extraction pattern mask
    for i=1:numel(struct)
         idx = strcmpi(struct{i}, EnDoFs{1});
         if any(idx)
             j = find(idx,1);
             if i == elementNo
                 %Degree of Freedom of the component to be extracted
                 nDoF = EnDoFs{2}{j};
                 pattern = [pattern, true(1, EnDoFs{2}{j})];
             else
                 pattern = [pattern, false(1, EnDoFs{2}{j})];
             end
         else
             fprintf('Element %s is not implemented!\n',struct{i});
             error('Unknown Element');
         end
    end
    %Get nb of triangles (the dirty way)
    nt = numel(vhIn)/numel(pattern);
    %Extraction mask
    mask = repmat(logical(pattern)',nt,1);
    vhMask = vhIn(mask);
    uOut = u(vhMask+1);
    vhOut = linspace(0,nDoF*nt-1,nDoF*nt);

end

function printhelp()
    fprintf('%s\n\n','Invalid call to ffreadfespace(). Correct usage is:');
    fprintf('%s\n',' -- [vhOut,uOut] = ffvectorget(struct,vhIn,u,elementNo)');
end
