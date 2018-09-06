classdef Mesh < handle
    % Class containing the coordinates and connectivities of the mesh
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
        mean_cell_size
        problem_characterisitc_length % !! Rename?? !!
    end
    
    methods
        function obj = create(obj,coordinates,connectivities)
            obj.coord = coordinates(:,1:obj.ndim);
            obj.connec = connectivities;
            obj.estimate_mesh_size;
            obj.estimate_mesh_characteristic_length;
        end
        
        function copy = duplicate(obj)
            copy = Mesh.create(obj.coord,obj.connec);
        end
    end
    
    methods (Access = private)
        function estimate_mesh_size(obj)
            x1 = obj.coord(obj.connec(:,1));
            x2 = obj.coord(obj.connec(:,2));
            x3 = obj.coord(obj.connec(:,3));
            
            x1x2 = abs(x2-x1);
            x2x3 = abs(x3-x2);
            x1x3 = abs(x1-x3);
            hs = max([x1x2,x2x3,x1x3]');
            
            obj.mean_cell_size = mean(hs);
        end
        
        function estimate_mesh_characteristic_length(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            obj.problem_characterisitc_length = norm(xmax-xmin)/2;
        end
    end
end

