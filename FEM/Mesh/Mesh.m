classdef Mesh < handle
    % Class containing the coordinates and connectivities of the mesh
    properties (GetAccess = public, SetAccess = protected)
        coord
        connec
        ndim
        geometryType
        mean_cell_size
        problem_characterisitc_length % !! Rename?? !!
    end
    
    methods
        function obj = create(obj,coordinates,connectivities)
            obj.coord = coordinates;
            obj.connec = connectivities;
            obj.ndim = size(coordinates,2);
            obj.setGeometryType;
            obj.estimate_mesh_size;
            obj.estimate_mesh_characteristic_length;
        end
        
        function copy = duplicate(obj)
            copy = Mesh;
            copy.create(obj.coord,obj.connec);
        end
    end
    
    methods (Access = private)
        function setGeometryType(obj)
            switch obj.ndim
                case 2
                    obj.setGeomType2D;
                case 3
                    obj.setGeomType3D;
            end
        end
        
        function setGeomType2D(obj)
            switch size(obj.connec,2)
                case 2
                    obj.geometryType = 'LINE';
                case 3
                    obj.geometryType = 'TRIANGLE';
                case 4
                    obj.geometryType = 'QUAD';
            end
        end
        
        function setGeomType3D(obj)
            switch size(obj.connec,2)
                case 4
                    obj.geometryType = 'TETRAHEDRA';
                case 8
                    obj.geometryType = 'HEXAHEDRA';
            end
        end
        
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

