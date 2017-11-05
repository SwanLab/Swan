classdef Mesh
    % Class containing the coordinates and connectivities of the mesh
    properties (GetAccess = public,SetAccess = private)
        % !! More elegant if Physical_Problem & subclasses !!
        nelem
        npnod
    end
    properties (GetAccess = {?Physical_Problem,?Geometry,?Postprocess,?TopOpt_Problem},SetAccess = private)
        % !! More elegant if Physical_Problem & subclasses !!
        coord
        connec
        geometryType
        ptype
        pdim
    end    
    
    methods
        function obj = Mesh(filename)
             
            data = Preprocess.readFromGiD(filename);
            obj.coord=data.xpoints(:,2:4);
            obj.connec=data.connectivities(:,2:length(data.connectivities(1,:)));
            obj.geometryType=data.geometry;
            obj.pdim=data.problem_dim;        
            obj.ptype=data.problem_type;
            obj.nelem = length(obj.connec(:,1));
            obj.npnod = length(obj.coord(:,1));
        end
    end
    
end

