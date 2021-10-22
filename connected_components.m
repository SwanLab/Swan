%% Local functions are part of the gptoolbox by Alec Jacobson
function C = connected_components(F)
% CONNECTED_COMPONENTS Determine the connected components of a mesh
% described by the simplex list F. Components are determined with respect
% to the edges of the mesh. That is, a single component may contain
% non-manifold edges and vertices.
%
% C = connected_components(F)
%
% Inputs:
%   F  #F by simplex-size list of simplices
% Outputs:
%   C  #V list of ids for each CC
%
% Examples:
%  trisurf(F,V(:,1),V(:,2),V(:,3), ...
%    connected_components([F;repmat(size(V,1),1,3)]));

% build adjacency list
A = adjacency_matrix(F);
[~,C] = conncomp(A);
end

function [A] = adjacency_matrix(E)
% ADJACENCY_MATRIX Build sparse adjacency matrix from edge list or face list
%
% [A] = adjacency_matrix(E)
% [A] = adjacency_matrix(F)
% [A] = adjacency_matrix(T)
%
% Inputs:
%   E  #E by 2 edges list
%   or
%   F  #F by 3 triangle list
%   or
%   T  #F by 4 tet list
% Outputs:
%   A  #V by #V adjacency matrix (#V = max(E(:)))
%
% See also: facet_adjacency_matrix
%

if size(E,2)>2
    F = E;
    E = meshEdges(F);
end

A = sparse([E(:,1) E(:,2)],[E(:,2) E(:,1)],1);
end

function [S,C] = conncomp(G)
% CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
% toobox. G is an n by n adjacency matrix, then this identifies the S
% connected components C. This is also an order of magnitude faster.
%
% [S,C] = conncomp(G)
%
% Inputs:
%   G  n by n adjacency matrix
% Outputs:
%   S  scalar number of connected components
%   C

% Transpose to match graphconncomp
G = G';

[p,~,r] = dmperm(G+speye(size(G)));
S = numel(r)-1;
C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
C(p) = C;
end

function edges = meshEdges(faces, varargin)
%MESHEDGES Computes array of edge vertex indices from face array
%
%   EDGES = meshEdges(FACES);
%
%   Example
%     meshEdges
%
%   See also
%     meshes3d, meshEdgeFaces, meshFaceEdges

% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-06-28,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

%   HISTORY
%   2013-08-22 rename from computeMeshEdges to meshEdges, add more control
%       on inputs

%% Process input arguments

if isstruct(faces) && isfield(faces, 'faces')
    % if input is a mesh structure, extract the 'faces' field
    faces = faces.faces;
elseif nargin > 1
    % if two arguments are given, keep the second one
    faces = varargin{1};
end


if ~iscell(faces)
    %% Process faces given as numeric array
    % all faces have same number of vertices, stored in nVF variable
    
    % compute total number of edges
    nFaces  = size(faces, 1);
    nVF     = size(faces, 2);
    nEdges  = nFaces * nVF;
    
    % create all edges (with double ones)
    edges = zeros(nEdges, 2);
    for i = 1:nFaces
        f = faces(i, :);
        edges(((i-1)*nVF+1):i*nVF, :) = [f' f([2:end 1])'];
    end
    
else
    %% faces are given as a cell array
    % faces may have different number of vertices
    
    % number of faces
    nFaces  = length(faces);
    
    % compute the number of edges
    nEdges = 0;
    for i = nFaces
        nEdges = nEdges + length(faces{i});
    end
    
    % allocate memory
    edges = zeros(nEdges, 2);
    ind = 0;
    
    % fillup edge array
    for i = 1:nFaces
        % get vertex indices, ensuring horizontal array
        f = faces{i}(:)';
        nVF = length(f);
        edges(ind+1:ind+nVF, :) = [f' f([2:end 1])'];
        ind = ind + nVF;
    end
    
end

% keep only unique edges, and return sorted result
edges = sortrows(unique(sort(edges, 2), 'rows'));

end


function varargout = removeMeshVertices(vertices, faces, indsToRemove)
%REMOVEMESHVERTICES Remove vertices and associated faces from a mesh
%
%   [V2, F2] = removeMeshVertices(VERTS, FACES, VERTINDS)
%   Removes the vertices specified by the vertex indices VERTINDS, and
%   remove the faces containing one of the removed vertices.
%
%   Example
%     % remove some vertices from a soccerball polyhedron
%     [v, f] = createSoccerBall; 
%     plane = createPlane([.6 0 0], [1 0 0]);
%     indAbove = find(~isBelowPlane(v, plane));
%     [v2, f2] = removeMeshVertices(v, f, indAbove);
%     drawMesh(v2, f2);
%     axis equal; hold on;
%
%   See also
%     meshes3d, trimMesh
 
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2016-02-03,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - Cepia Software Platform.

% parse inputs
if nargin == 2
    indsToRemove = faces;
    [vertices, faces] = parseMeshData(vertices);
end

% create array of indices to keep
nVertices = size(vertices, 1);
newInds = (1:nVertices)';
newInds(indsToRemove) = [];

% create new vertex array
vertices2 = vertices(newInds, :);

% compute map from old indices to new indices
oldNewMap = zeros(nVertices, 1);
for iIndex = 1:size(newInds, 1)
   oldNewMap(newInds(iIndex)) = iIndex; 
end

% change labels of vertices referenced by faces
if isnumeric(faces)
    faces2 = oldNewMap(faces);
    if size(faces2,2)==1; faces2=faces2'; end
    % keep only faces with valid vertices
    faces2 = faces2(sum(faces2 == 0, 2) == 0, :);
elseif iscell(faces)
    faces2 = cell(1, length(faces));
    for iFace = 1:length(faces)
         newFace = oldNewMap(faces{iFace});
         % add the new face only if all vertices are valid
         if ~any(newFace == 0)
             faces2{iFace} = newFace;
         end
    end
    
    % remove empty faces
    faces2 = faces2(~cellfun(@isempty, faces2));
end

% format output arguments
varargout = formatMeshOutput(nargout, vertices2, faces2);

end

function res = formatMeshOutput(nbArgs, vertices, edges, faces)
%FORMATMESHOUTPUT Format mesh output depending on nargout
%
%   OUTPUT = formatMeshOutput(NARGOUT, VERTICES, EDGES, FACES)
%   Utilitary function to convert mesh data .
%   If NARGOUT is 0 or 1, return a matlab structure with fields vertices,
%   edges and faces.
%   If NARGOUT is 2, return a cell array with data VERTICES and FACES.
%   If NARGOUT is 3, return a cell array with data VERTICES, EDGES and
%   FACES. 
%
%   OUTPUT = formatMeshOutput(NARGOUT, VERTICES, FACES)
%   Same as before, but do not intialize EDGES in output. NARGOUT can not
%   be equal to 3.
%
%   Example
%     % Typical calling sequence (for a very basic mesh of only one face)
%     v = [0 0; 0 1;1 0;1 1];
%     e = [1 2;1 3;2 4;3 4];
%     f = [1 2 3 4];
% 
%     varargout = formatMeshOutput(nargout, v, e, f);
%
%   See also
%   meshes3d, parseMeshData
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2010-12-06,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2010 INRA - Cepia Software Platform.

if nargin < 4
    faces = edges;
    edges = [];
end

switch nbArgs
    case {0, 1}
        % output is a data structure with fields vertices, edges and faces
        mesh.vertices = vertices;
        mesh.edges = edges;
        mesh.faces = faces;
        res = {mesh};

    case 2
        % keep only vertices and faces
        res = cell(nbArgs, 1);
        res{1} = vertices;
        res{2} = faces;
        
    case 3
        % return vertices, edges and faces as 3 separate outputs
        res = cell(nbArgs, 1);
        res{1} = vertices;
        res{2} = edges;
        res{3} = faces;
        
    otherwise
        error('Can not manage more than 3 outputs');
end
end