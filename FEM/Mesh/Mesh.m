classdef Mesh < handle
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
        ndim
        geometryType
        %         mean_cell_size
        problem_characterisitc_length % !! Rename?? !!
    end
    
    methods (Access = public)
        function obj = create(obj,coordinates,connectivities)
            obj.coord = coordinates;
            obj.connec = connectivities;
            obj.ndim = size(coordinates,2);
            obj.computeGeometryType;
            %             obj.mean_cell_size = obj.computeMeanCellSize;
            obj.estimate_mesh_characteristic_length;
        end
        
        function copy = clone(obj)
            copy = Mesh;
            copy.create(obj.coord,obj.connec);
        end
        
        function meanCellSize = computeMeanCellSize(obj)
            x1 = obj.coord(obj.connec(:,1));
            x2 = obj.coord(obj.connec(:,2));
            x3 = obj.coord(obj.connec(:,3));
            
            x1x2 = abs(x2-x1);
            x2x3 = abs(x3-x2);
            x1x3 = abs(x1-x3);
            hs = max([x1x2,x2x3,x1x3]');
            
            meanCellSize = mean(hs);
        end
    end
    
    methods (Access = protected)
        function computeGeometryType(obj)
            nnode = size(obj.connec,2);
            obj.geometryType = MeshGeometryType_Factory.getGeometryType(obj.ndim,nnode);
        end
    end
    
    methods (Access = private)
        
        
        function estimate_mesh_characteristic_length(obj)
            xmin = min(obj.coord);
            xmax = max(obj.coord);
            obj.problem_characterisitc_length = norm(xmax-xmin)/2;
        end
    end
end

