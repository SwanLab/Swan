% Vectorized version of createPeriodicMesh2D_withBoundary
function [COOR_ALL, CN_ALL, CNb_ALL, MAT_ALL, ELEMENTS_PER_COPY] = ...
    createPeriodicMesh2D_withBoundary_CGPT_vect(COOR, CN, CNb, px, py, MAT)
%--------------------------------------------------------------------------
% Tiling of a 2D unit cell mesh including both interior and boundary elements.
% Vectorized implementation.
%--------------------------------------------------------------------------

% Parameters
tol = 1e-4;
makeKey = @(coord) sprintf('%d_%d', round(coord(1)/tol), round(coord(2)/tol));

% Dimensions
[Nnodes, dim] = size(COOR);
[Nelements, nnodeE] = size(CN);
[Nbndelems, nnodeB] = size(CNb);

% Base shifts
dx = max(COOR(:,1)) - min(COOR(:,1));
dy = max(COOR(:,2)) - min(COOR(:,2));

% Initialize containers
node_counter = 0;
coord_map = containers.Map('KeyType','char', 'ValueType','int32');
COOR_dict = zeros(px*py*Nnodes, 2);

CN_ALL = zeros(px*py*Nelements, nnodeE);
CNb_ALL = zeros(px*py*Nbndelems, nnodeB);
MAT_ALL = zeros(px*py*Nelements, 1);
ELEMENTS_PER_COPY = cell(px*py,1);

% Loop over tiles (ix, iy)
tile_id = 0;
elem_counter = 0;
bnd_counter = 0;

[X, Y] = meshgrid(0:px-1, 0:py-1);
shifts = [X(:)*dx, Y(:)*dy];

for t = 1:size(shifts, 1)
    tile_id = tile_id + 1;
    shift = shifts(t,:);
    COOR_shifted = COOR + shift;

    % Vectorized coordinate key mapping
    key_array = arrayfun(@(i) makeKey(COOR_shifted(i,:)), 1:Nnodes, 'UniformOutput', false);

    idx_new = zeros(Nnodes, 1);
    for inode = 1:Nnodes
        key = key_array{inode};
        if ~isKey(coord_map, key)
            node_counter = node_counter + 1;
            coord_map(key) = node_counter;
            COOR_dict(node_counter,:) = COOR_shifted(inode,:);
        end
        idx_new(inode) = coord_map(key);
    end

    % Build element connectivity
    CN_tile = idx_new(CN);
    CNb_tile = idx_new(CNb);

    elem_range = elem_counter + (1:Nelements);
    bnd_range = bnd_counter + (1:Nbndelems);

    CN_ALL(elem_range,:) = CN_tile;
    CNb_ALL(bnd_range,:) = CNb_tile;
    MAT_ALL(elem_range) = MAT;
    ELEMENTS_PER_COPY{tile_id} = elem_range(:);

    elem_counter = elem_counter + Nelements;
    bnd_counter = bnd_counter + Nbndelems;
end

COOR_ALL = COOR_dict(1:node_counter,:);
end
