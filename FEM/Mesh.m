classdef Mesh
    % Class containing the coordinates and connectivities of the mesh
    properties (GetAccess = public,SetAccess = private)
        % !! More elegant if Physical_Problem & subclasses !!
        nelem
        npnod
        coord
        connec
        pdim
    end
    
    properties (GetAccess = public,SetAccess = public)
        % !! More elegant if Physical_Problem & subclasses !!
        ptype
        scale
    end
    
    
    properties (GetAccess = {?Physical_Problem,?Geometry,?Postprocess,?TopOpt_Problem},SetAccess = {?Physical_Problem, ?Element_DiffReact}) % !! Element_DiffReact -> Chapusilla !!
        % !! More elegant if Physical_Problem & subclasses !!
        geometryType

    end
    
    
    methods (Access = ?Physical_Problem)
        function obj = Mesh(filename)
            data = Preprocess.readFromGiD(filename);
            obj.coord = data.xpoints(:,2:4);
            obj.connec = data.connectivities(:,2:length(data.connectivities(1,:)));
            obj.geometryType = data.geometry;
            obj.pdim = data.problem_dim;
            obj.ptype = data.problem_type;
            obj.nelem = length(obj.connec(:,1));
            obj.npnod = length(obj.coord(:,1));
            obj.scale = data.scale;
        end
    end
    
end

