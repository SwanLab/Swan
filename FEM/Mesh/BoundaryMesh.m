classdef BoundaryMesh < handle
    
    properties (Access = public)
       nodesInBoxFaces
       globalConnec
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
            obj.createGlobalConnec();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nodesInBoxFaces = cParams.nodesInBoxFaces;
            obj.connec          = cParams.connec;
            obj.coord           = cParams.coord;
        end
        
        function createMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            s.kFace  = -1;
            obj.mesh = Mesh(s);
        end
        
        function createGlobalConnec(obj)
            nodes = find(obj.nodesInBoxFaces);
            for inode = 1:size(obj.connec,2)
               obj.globalConnec(:,inode) = nodes(obj.connec(:,inode));
            end            
        end
    end
    
end