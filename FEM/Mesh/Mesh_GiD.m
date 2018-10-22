classdef Mesh_GiD < Mesh
    % Class containing the coordinates and connectivities of the mesh
    properties (GetAccess = public,SetAccess = public)
        pdim
        dirichlet
        pointload
    end
    
    properties (GetAccess = public,SetAccess = public)
        ptype
        scale
        problemID
    end
    
    methods
        function obj = Mesh_GiD(filename)
            if nargin > 0
                data = Preprocess.readFromGiD(filename);
                obj.pdim = data.problem_dim;
                obj.getNdim;                
                obj.geometryType = data.geometry;
                obj.ptype = data.problem_type;
                obj.scale = data.scale;
                obj.problemID=filename;
                obj.create(data.xpoints(:,2:obj.ndim+1),data.connectivities(:,2:end));
                
                if strcmpi(data.problem_type,'elastic')
                    obj.dirichlet = data.dirichlet_data;
                    obj.pointload = data.pointload;
                end
            end
        end
        
        function copy = duplicate(obj)
            copy = Mesh_GiD(obj.problemID);
        end
        
        function simplified_copy = getSimplifiedMesh(obj)
            simplified_copy = Mesh;
            simplified_copy.create(obj.coord,obj.connec);
        end
    end
    
    methods (Access = private)
        function getNdim(obj)
            switch obj.pdim
                case '2D'
                    obj.ndim = 2;
                case '3D'
                    obj.ndim = 3;
            end
        end
    end
end

