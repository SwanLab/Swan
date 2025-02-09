function mshWriteMsh(varargin)
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
%|    #    |   FILE       : mshWriteMsh.m                                 |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF : 31.12.2018                                    |
%| ( === ) |   SYNOPSIS   : Write mesh and data to msh format             |
%|  `---'  |                (particle, edge, triangular and tetrahedral)  |
%+========================================================================+

% Input analysis
filename = varargin{1};
mesh     = varargin{2};
if (nargin == 3)
    data = varargin{3};
end

% Security
if (size(mesh.vtx,2) ~= 3) || (size(mesh.elt,2) > 4)
    error('mshWriteMsh.m : unavailable case')
end
   
% Open file
fid = fopen(filename,'w');

% Header
fprintf(fid,'%s\n','$MeshFormat');
fprintf(fid,'%s\n','2.2 0 8');
fprintf(fid,'%s\n','$EndMeshFormat');

% Nodes
fprintf(fid,'%s\n','$Nodes');
fprintf(fid,'%d\n',size(mesh.vtx,1));
for i = 1:size(mesh.vtx,1)
    fprintf(fid,'%d %f %f %f\n',i,mesh.vtx(i,:));
end
fprintf(fid,'%s\n','$EndNodes');

% Elements (up to tetra)
fprintf(fid,'%s\n','$Elements');
fprintf(fid,'%d\n',size(mesh.elt,1));
for i = 1:size(mesh.elt,1)
    if (size(mesh.elt,2) == 1)  % particle
        fprintf(fid,'%d %d %d %d %d  %d\n',i,15,2,mesh.col(i),0,mesh.elt(i,:));
    elseif (size(mesh.elt,2) == 2)  % segment
        fprintf(fid,'%d %d %d %d %d  %d %d\n',i,1,2,mesh.col(i),0,mesh.elt(i,:));
    elseif (size(mesh.elt,2) == 3)  % triangle
        fprintf(fid,'%d %d %d %d %d  %d %d %d\n',i,2,2,mesh.col(i),0,mesh.elt(i,:));
    elseif (size(mesh.elt,2) == 4)  % tetra
        fprintf(fid,'%d %d %d %d %d  %d %d %d %d\n',i,4,2,mesh.col(i),0,mesh.elt(i,:));
    end
end
fprintf(fid,'%s\n','$EndElements');

% Nodes data
if (nargin == 3)
    n = length(data(1,:));        % 1:scalar, 3:vector, 9:tensor
    fprintf(fid,'%s\n','$NodeData');
    fprintf(fid,'%d\n',1);        % number-of-string-tag # must be 1 for Mmg
    fprintf(fid,'%s:metric\n',filename); % string # name of the metric field
    fprintf(fid,'%d\n',1);        % number-of-real-tags # must be 1 for Mmg
    fprintf(fid,'%d\n',0);        % real # ignored value
    fprintf(fid,'%d\n',3);        % number-of-integer-tags # must be 3 for Mmg  
    fprintf(fid,'%d\n',0);        % ignored value
    fprintf(fid,'%d\n',n);        % metric-type # type of metric  
    fprintf(fid,'%d\n',size(mesh.vtx,1)); % number of metrics: must match with the number of nodes in the $Nodes field
    for i = 1:size(mesh.vtx,1)
        if (n == 1)
            fprintf(fid,'%d %f\n',i,data(i,:));
        elseif (n == 3)
            fprintf(fid,'%d %f %f %f\n',i,data(i,:));
        elseif (n == 9)
            fprintf(fid,'%d %f %f %f %f %f %f %f %f %f\n',i,data(i,:));
        else
            error('mshWriteMsh.m : unavailable case');
        end
    end
    fprintf(fid,'%s\n','$EndNodeData');
end

% Close file
fclose(fid);
disp([filename,' created.']);
end
