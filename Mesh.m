classdef Mesh
    % Class containing the coordinates and connectivities of the mesh
    
    properties
        coord
        connec
        nelem
        npnod
    end
    
    
    methods
        function obj = Mesh()
            [obj.coord,obj.connec] = Preprocess.readFromGiD();
            obj.coord = obj.coord(:,2:4);
            obj.connec = obj.connec(:,2:4);
            obj.nelem = length(obj.connec(:,1));
            obj.npnod = length(obj.coord(:,1));
        end
    end
    
end

