function mshWriteVtk(varargin)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : mshWriteVtk.m                                 |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Write mesh and data to vtk format             |
%|  `---'  |                (only surfacic triangular mesh)               |
%+========================================================================+

% Input anamysis
filename = varargin{1};
mesh     = varargin{2};
if (nargin == 2)
    V = [];
elseif (nargin == 3)
    V = varargin{3};
end

% Security
if (size(mesh.vtx,2) ~= 3) || (size(mesh.elt,2) ~= 3)
    error('mshWriteVtk.m : unavailable case')
end
   
% Open file
fid = fopen(filename,'w');

% Header
fprintf(fid,'%s\n','# vtk DataFile Version 3.0 #');
fprintf(fid,'%s\n','# COMPUTE WITH Gypsilab 0.32 #');
fprintf(fid,'%s\n','ASCII');
fprintf(fid,'%s\n','DATASET UNSTRUCTURED_GRID');

% Node list
fprintf(fid,'%s %d %s\n','POINTS',size(mesh.vtx,1),'FLOAT');
for i = 1:size(mesh.vtx,1)
    fprintf(fid,'%f %f %f\n',mesh.vtx(i,:));
end

% Element list
fprintf(fid,'%s %d %d\n','CELLS',size(mesh.elt,1),4*size(mesh.elt,1));
for i = 1:size(mesh.elt,1)
    fprintf(fid,'%d  %d  %d  %d\n',3,mesh.elt(i,:)-1);
end

% Element type triangle
type_VTK = 5;
fprintf(fid,'%s %d\n','CELL_TYPES',size(mesh.elt,1));
for i = 1:size(mesh.elt,1)
    fprintf(fid,'%d\n',type_VTK);
end

% Vertex data
if ~isempty(V) && (size(V,1) == size(mesh.vtx,1))
    fprintf(fid,'%s %d\n','POINT_DATA',size(mesh.vtx,1));
    fprintf(fid,'%s \n', 'SCALARS (REAL,IMAG,0) FLOAT 3');
    fprintf(fid,'%s \n', 'LOOKUP_TABLE default');
    for i = 1:size(mesh.vtx,1)
        fprintf(fid,'%f %f %f\n',real(V(i)),imag(V(i)),0);
    end
end

% Element data
if ~isempty(V) && (size(V,1) == size(mesh.elt,1))
    fprintf(fid,'%s %d\n','CELL_DATA',size(mesh.elt,1));
    fprintf(fid,'%s \n','SCALARS (REAL,IMAG,0) FLOAT 3');
    fprintf(fid,'%s \n','LOOKUP_TABLE default');
    for i = 1:size(mesh.elt,1)
        fprintf(fid,'%f %f %f\n',real(V(i)),imag(V(i)),0);
    end
end

% Close file
fclose(fid);
disp([filename,' created.']);
end
