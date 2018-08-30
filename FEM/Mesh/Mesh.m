classdef Mesh
    % Class containing the coordinates and connectivities of the mesh
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
    end

    methods
        function obj = create(coordinates,connectivities)
            obj.coord = coordinates;
            obj.connec = connectivities;
        end
    end
end

