
function [COOR_ALL, CN_ALL, CNb_ALL, MAT_ALL, ELEMENTS_PER_COPY] = ...
    createPeriodicMesh2D_withBoundary(COOR, CN, CNb, px, py, MAT)
%--------------------------------------------------------------------------
% Tiling of a 2D unit cell mesh including both interior and boundary elements.
%
% INPUT:
%   COOR : [Nnodes x 2]         - coordinates of unit cell
%   CN   : [Nelements x nnodeE] - connectivity of interior elements
%   CNb  : [Nbndelems x nnodeB] - connectivity of boundary elements
%   px   : number of repetitions in x-direction
%   py   : number of repetitions in y-direction
%   MAT  : [Nelements x 1]      - material index per element
%
% OUTPUT:
%   COOR_ALL          : unique global coordinates
%   CN_ALL            : global interior connectivity
%   CNb_ALL           : global boundary connectivity
%   MAT_ALL           : global material indices
%   ELEMENTS_PER_COPY : cell array of global element indices per tile
%--------------------------------------------------------------------------

tol = 1e-4;
makeKey = @(coord) sprintf('%d_%d', round(coord(1)/tol), round(coord(2)/tol));

% Dimensions
Nnodes    = size(COOR, 1);
Nelements = size(CN, 1);
Nbndelems = size(CNb, 1);
nnodeE    = size(CN, 2);
nnodeB    = size(CNb, 2);

dx = max(COOR(:,1)) - min(COOR(:,1));
dy = max(COOR(:,2)) - min(COOR(:,2));

% Storage
node_counter = 0;
coord_map = containers.Map('KeyType','char', 'ValueType','int32');
COOR_dict = [];

elem_counter = 0;
bnd_counter = 0;
CN_ALL = [];
CNb_ALL = [];
MAT_ALL = [];
ELEMENTS_PER_COPY = cell(px*py,1);

tile_id = 0;
disp('Creating tiled mesh')
for iy = 0:py-1
    
    for ix = 0:px-1
        disp(['Domain ',num2str(iy+1),' ',num2str(ix+1)])
        tile_id = tile_id + 1;
        shift = [ix*dx, iy*dy];
        COOR_shifted = COOR + shift;
        
        CN_tile = zeros(Nelements, nnodeE);
        CNb_tile = zeros(Nbndelems, nnodeB);
        
        for inode = 1:Nnodes
            coord = COOR_shifted(inode,:);
         %   key = sprintf('%.10f_%.10f', coord(1), coord(2));
            key = makeKey(coord); 
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
              %  coord = round(coord / tol) * tol;  % round to nearest tol
               % key = sprintf('%.10f_%.10f', coord(1), coord(2));
                 key = makeKey(coord); 
                
                %  key = sprintf('%.10f_%.10f', coord(1), coord(2));
                CN_tile(ie,j) = coord_map(key);
            end
        end
        
        for ib = 1:Nbndelems
            for j = 1:nnodeB
                local_node = CNb(ib,j);
                coord = COOR_shifted(local_node,:);
                
              %  coord = round(coord / tol) * tol;  % round to nearest tol
               key = makeKey(coord); 
              %  key = sprintf('%.10f_%.10f', coord(1), coord(2));
                
                %   key = sprintf('%.10f_%.10f', coord(1), coord(2));
                CNb_tile(ib,j) = coord_map(key);
            end
        end
        
        elem_range = elem_counter + (1:Nelements);
        bnd_range = bnd_counter + (1:Nbndelems);
        
        CN_ALL(elem_range,:) = CN_tile;
        CNb_ALL(bnd_range,:) = CNb_tile;
        MAT_ALL(elem_range,1) = MAT;
        ELEMENTS_PER_COPY{tile_id} = elem_range(:);
        
        elem_counter = elem_counter + Nelements;
        bnd_counter = bnd_counter + Nbndelems;
    end
end

COOR_ALL = COOR_dict;

end