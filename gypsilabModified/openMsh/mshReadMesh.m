function [vtx,elt,col] = mshReadMesh(filename)
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
%|    #    |   FILE       : mshReadMesh.m                                 |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Read .mesh files (vertex, triangular elements,|
%|  `---'  |                and colours)                                  |
%+========================================================================+

% Open file
fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
end

% Nodes
str = fgets(fid);
while isempty(strfind(str,'Vertices'))
    str = fgets(fid);
end
Nvtx = str2double(fgets(fid));
vtx  = zeros(Nvtx,3);
for i = 1:Nvtx
    tmp = str2num(fgets(fid));
    vtx(i,:) = tmp(1:3);
end

% Support only triangular elements
str = fgets(fid);
if strcmp(str(1:9),'Triangles') 
    Nelt = str2double(fgets(fid));
    elt  = zeros(Nelt,4);
    for i = 1:Nelt
        elt(i,:) = str2num(fgets(fid));
    end
    col = elt(:,4);
    elt = elt(:,1:3); 
    
elseif strcmp(str(1:10),'Tetrahedra')
    Nelt = str2double(fgets(fid));
    elt  = zeros(Nelt,5);
    for i = 1:Nelt
        elt(i,:) = str2num(fgets(fid));
    end
    col = elt(:,5);
    elt = elt(:,1:4);
    
else
    error('mshReadMesh : unavailable case');
end

% Close file
fclose(fid);
end