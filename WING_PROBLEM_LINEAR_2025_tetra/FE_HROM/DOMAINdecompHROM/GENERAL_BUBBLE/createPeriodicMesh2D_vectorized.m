function [COOR_ALL, CN_ALL, MAT_ALL, ELEMENTS_PER_COPY] = ...
    createPeriodicMesh2D_vectorized(COOR, CN, px, py, MAT)
%--------------------------------------------------------------------------
% Tiling of a 2D unit cell mesh with full preservation of:
% - element order
% - node order in each element
% - node IDs remapped only due to deduplication (not reordered)
%
% INPUT:
%   COOR : [Nnodes x 2]         - coordinates of unit cell
%   CN   : [Nelements x nnodeE] - connectivity of unit cell
%   px   : number of repetitions in x-direction
%   py   : number of repetitions in y-direction
%   MAT  : [Nelements x 1]      - material index per element
%
% OUTPUT:
%   COOR_ALL          : unique global coordinates
%   CN_ALL            : global connectivity (order preserved)
%   MAT_ALL           : global material indices
%   ELEMENTS_PER_COPY : cell array of global element indices per tile
%--------------------------------------------------------------------------

tol = 1e-10;

% Dimensions
Nnodes    = size(COOR, 1);
Nelements = size(CN, 1);
nnodeE    = size(CN, 2);

dx = max(COOR(:,1)) - min(COOR(:,1));
dy = max(COOR(:,2)) - min(COOR(:,2));

% Storage
COOR_all = [];
CN_all = [];
MAT_ALL = [];
ELEMENTS_PER_COPY = cell(px*py,1);

node_counter = 0;
coord_map = containers.Map('KeyType','char', 'ValueType','int32'); % coordinate → global index
COOR_dict = [];

elem_counter = 0;
CN_ALL = [];

tile_id = 0;
for iy = 0:py-1
    for ix = 0:px-1
        tile_id = tile_id + 1;
        shift = [ix*dx, iy*dy];
        COOR_shifted = COOR + shift;

        CN_tile = zeros(Nelements, nnodeE);

        for inode = 1:Nnodes
            coord = COOR_shifted(inode,:);
            key = sprintf('%.10f_%.10f', coord(1), coord(2));
            if ~isKey(coord_map, key)
                node_counter = node_counter + 1;
                coord_map(key) = node_counter;
                COOR_dict(node_counter,:) = coord;
            end
        end

        for ie = 1:Nelements
            for j = 1:nnodeE
                local_node = CN(ie,j);
                coord = COOR_shifted(local_node,:);
                key = sprintf('%.10f_%.10f', coord(1), coord(2));
                CN_tile(ie,j) = coord_map(key);
            end
        end

        elem_range = elem_counter + (1:Nelements);
        CN_ALL(elem_range,:) = CN_tile;
        MAT_ALL(elem_range,1) = MAT;
        ELEMENTS_PER_COPY{tile_id} = elem_range(:);
        elem_counter = elem_counter + Nelements;
    end
end

COOR_ALL = COOR_dict;

end
