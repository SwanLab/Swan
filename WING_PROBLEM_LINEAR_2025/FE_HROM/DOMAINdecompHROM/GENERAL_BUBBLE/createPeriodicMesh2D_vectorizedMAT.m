function [COOR_ALL, CN_ALL, MAT_ALL, ELEMENTS_PER_COPY] = ...
    createPeriodicMesh2D_vectorizedMAT(COOR, CN, px, py, MAT)
%--------------------------------------------------------------------------
% Create a periodic 2D mesh by tiling a unit cell mesh.
% Merges coinciding nodes, assigns material indices, and maps elements to each copy.
%
% INPUT:
%   COOR : Nnodes × 2 → Coordinates of unit cell
%   CN   : Nelements × nnodeE → Connectivity of unit cell
%   px   : Repetitions in x direction
%   py   : Repetitions in y direction
%   MAT  : Nelements × 1 → Material index per unit cell element
%
% OUTPUT:
%   COOR_ALL          : Global coordinates after node merging
%   CN_ALL            : Global connectivity
%   MAT_ALL           : Global material index vector
%   ELEMENTS_PER_COPY : Cell array with element indices per tile copy
%  JAHO, 28-May-2025, guided by ChatGPT4
%--------------------------------------------------------------------------
% Compute unit cell size
xL = max(COOR(:,1)) - min(COOR(:,1));
yL = max(COOR(:,2)) - min(COOR(:,2));

% Sizes
Nnodes = size(COOR,1);
Nelements = size(CN,1);
nnodeE = size(CN,2);

% All shifts (px * py)
[Ix, Iy] = meshgrid(0:px-1, 0:py-1);
Ix = Ix(:); Iy = Iy(:);
nCopies = numel(Ix);

% Initialize
COOR_all = zeros(Nnodes * nCopies, 2);
CN_all = zeros(Nelements * nCopies, nnodeE);
MAT_ALL = zeros(Nelements * nCopies, 1);
ELEMENTS_PER_COPY = cell(nCopies,1);

for k = 1:nCopies
    % Node translation
    delta = [Ix(k)*xL, Iy(k)*yL];
    node_range = (k-1)*Nnodes + (1:Nnodes);
    elem_range = (k-1)*Nelements + (1:Nelements);
    
    COOR_all(node_range, :) = COOR + delta;
    CN_all(elem_range, :) = CN + (k-1)*Nnodes;
    MAT_ALL(elem_range) = MAT;
    ELEMENTS_PER_COPY{k} = elem_range(:);
end

% Collapse duplicate nodes
[COOR_ALL, ~, ic] = unique(round(COOR_all, 10), 'rows');

CN_ALL = ic(CN_all);  % No reshape needed, shape and order preserved

% 
% CN_ALL = reshape(ic(CN_all), [], nnodeE);
end
