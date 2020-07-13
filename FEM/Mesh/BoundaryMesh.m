classdef BoundaryMesh < handle
    
    properties (Access = public)
       nodesInBoxFaces
       globalConnectivities
       mesh
    end
    
    properties (Access = private)
       connec
       coord
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = BoundaryMesh(cParams)
            obj.init(cParams)
            obj.createMesh();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodesInBoxFaces = cParams.nodesInBoxFaces;
            obj.globalConnec    = cParams.globalConnec;
            obj.connec          = cParams.connec;
            obj.coord           = cParams.coord;
        end
        
        function createMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            s.type   = 'BOUNDARY';
            obj.mesh = Mesh().create(s);
        end
    end
    
end