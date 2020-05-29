function mshWritePly(filename,mesh)
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
%|    #    |   FILE       : mshWritePly.m                                 |
%|    #    |   VERSION    : 0.43                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2018                                    |
%| ( === ) |   SYNOPSIS   : Write mesh and data to ply format             |
%|  `---'  |                (only surfacic triangular mesh)               |
%+========================================================================+

% Security
if (size(mesh.vtx,2) ~= 3) || (size(mesh.elt,2) ~= 3)
    error('mshWritePly.m : unavailable case')
end
   
% Open file
fid = fopen(filename,'w');

% Header
fprintf(fid,'%s\n','ply');
fprintf(fid,'%s\n','format ascii 1.0');
fprintf(fid,'%s\n','comment Gypsilab generated PLY file');
fprintf(fid,'%s %d\n','element vertex ',size(mesh.vtx,1));
fprintf(fid,'%s\n','property float x');
fprintf(fid,'%s\n','property float y');
fprintf(fid,'%s\n','property float z');
fprintf(fid,'%s %d\n','element face ',size(mesh.elt,1));
fprintf(fid,'%s\n','property list uchar int vertex_indices');
fprintf(fid,'%s\n','end_header');

% Node list
for i = 1:size(mesh.vtx,1)
    fprintf(fid,'%f %f %f\n',mesh.vtx(i,:));
end

% Element list
for i = 1:size(mesh.elt,1)
    fprintf(fid,'%d  %d  %d  %d\n',3,mesh.elt(i,:)-1);
end

% Close file
fclose(fid);
disp([filename,' created.']);
end
