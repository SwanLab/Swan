function mmgMeshWrite(varargin)
% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique
% Algiane Froehly, CARDAMOME, INRIA-SOFT
% LGPL Lesser General Public License v3.0.
% Remeshing using Mmg tools : https://www.mmgtools.org

% Input analysis
filename = varargin{1};
mesh     = varargin{2};
if (nargin==3)
    I = varargin{3};    
elseif (nargin==5)
    I = varargin{3};
    RequieredEdges = varargin{4};
    edges = varargin{5};
else 
    I = [];
end

% Space dimension 
if is2d(mesh)
    dim = 2;
else
    dim = 3;
end

% Open .mesh file
fid = fopen(filename,'w');
if(fid==-1)
    error('mmgMeshWrite.m : unavailable case');
end

% Header
fprintf(fid,'%s\n','# MESH generated using Gyspilab for Matlab #');
fprintf(fid,'%s\n\n','MeshVersionFormatted 2');
fprintf(fid,'%s\n\n',['Dimension ',num2str(dim)]);

% Write nodes
fprintf(fid,'Vertices\n');
N = size(mesh.vtx,1);
fprintf(fid,'%d\n',N);
if is2d(mesh)
    fprintf(fid,'%f %f %d\n',[mesh.vtx(:,1:2),zeros(N,1)]');
else
    fprintf(fid,'%f %f %f %d\n',[mesh.vtx,zeros(N,1)]');
end
fprintf(fid,'\n');

% Write elements
if (size(mesh.elt,2)==3)
    fprintf(fid,'Triangles\n');
elseif (size(mesh.elt,2)==4)
    fprintf(fid,'Tetrahedra\n');
else
    error('mmgMeshWrite.m : unavailable case');
end
fprintf(fid,'%d\n',size(mesh.elt,1));
fprintf(fid,'%d %d %d %d\n',[mesh.elt,mesh.col]');
fprintf(fid,'\n');

% Write required elements
if ~isempty(I) 
    if (size(mesh.elt,2)==3)
        fprintf(fid,'RequiredTriangles\n');
    elseif (size(mesh.elt,2)==4)
        fprintf(fid,'RequiredTetrahedra\n');
    end
    fprintf(fid,'%d\n',length(I));
    fprintf(fid,'%d\n',I');
    fprintf(fid,'\n');
end

% msRelator = MasterSlaveRelator(mesh.vtx(:,1:2));
% master_slave = msRelator.getRelation();
% 
% 



if (nargin==5)
    
% Write nodes
fprintf(fid,'Edges\n');
N = size(edges,1);
fprintf(fid,'%d\n',N);
if is2d(mesh)
    fprintf(fid,'%d %d %d\n',edges');
else
    fprintf(fid,'%f %f %f %d\n',[mesh.vtx,zeros(N,1)]');
end
fprintf(fid,'\n');    

fprintf(fid,'RequiredEdges\n');
fprintf(fid,'%d\n',length(RequieredEdges));
fprintf(fid,'%d\n',RequieredEdges);
fprintf(fid,'\n');
end
% Termination
fprintf(fid,'End \n');

% Close mesh file
fclose(fid);
end
