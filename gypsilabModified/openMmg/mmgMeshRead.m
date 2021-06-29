function [mesh,req,label,edg] = mmgMeshRead(filename)
% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique
% Algiane Froehly, CARDAMOME, INRIA-SOFT
% LGPL Lesser General Public License v3.0.
% Remeshing using Mmg tools : https://www.mmgtools.org

% Open .mesh file
fid = fopen(filename,'r');
if (fid==-1)
    error('mmgMeshRead.m : unavailable case');
end

% Initialization
dim  = [];
vtx  = [];
tri  = [];
tet  = [];
Itri = [];
Itet = [];

% Read file using keywords
str = fgets(fid);
while ~(str==-1) 
    % Dimension
    if startsWith(str,'Dimension')
        dim = str2double(str(end-1));
    
    % Vertices
    elseif startsWith(str,'Vertices')
        N = str2double(fgets(fid));
        if (dim==2)
            vtx = fscanf(fid,'%f %f %d',[3 N])';
        elseif (dim==3)
            vtx = fscanf(fid,'%f %f %f %d',[4 N])';
        end
        
    % Triangular elements and colours
    elseif startsWith(str,'Triangles')
        N   = str2double(fgets(fid));
        tri = fscanf(fid,'%d %d %d %d',[4 N])';
    
    % Tetrahedral elements and colours
    elseif startsWith(str,'Tetrahedra')
        N   = str2double(fgets(fid));
        tet = fscanf(fid,'%d %d %d %d %d',[5 N])';
        
    % Required triangles
    elseif startsWith(str,'RequiredTriangles')
        N    = str2double(fgets(fid));
        Itri = fscanf(fid,'%d',[1 N])';
        
    % Required tetrahedra
    elseif startsWith(str,'RequiredTetrahedra')
        N    = str2double(fgets(fid));
        Itet = fscanf(fid,'%d',[1 N])';
   
    % Triangular elements and colours
    elseif startsWith(str,'Edges')
        N   = str2double(fgets(fid));
        edg = fscanf(fid,'%d %d %d',[3 N])';   
    end
    
    % Next line
    str = fgets(fid);
end

% Close file
fclose(fid);

% Final vertices 
label = [];
if (dim==2) 
    label = edg(:,3);
    vtx(:,3) = 0;
end
vtx = vtx(:,1:3);

% Final elements (tetra > triangle > to do)
if ~isempty(tet)
    elt = tet(:,1:4);
    col = tet(:,5);
    req = Itet;
elseif ~isempty(tri)
    elt = tri(:,1:3);
    col = tri(:,4);
    req = Itri;
else
    error('mmgMeshRead.m : unavailable case');
end

% Final mesh
mesh = msh(vtx,elt,col);
end
