classdef BoundaryMesh < handle
    
    properties (Access = public)
       nodesInBoxFaces
       mesh
       globalConnec
       dimension
       isRectangularBox
    end
    
    properties (Access = private)
       connec
       coord
       kFace
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
            obj.nodesInBoxFaces  = cParams.nodesInBoxFaces;
            obj.connec           = cParams.connec;
            obj.coord            = cParams.coord;
            obj.dimension        = cParams.dimension;
            obj.kFace            = cParams.kFace;
            obj.isRectangularBox = cParams.isRectangularBox;
        end
        
        function createMesh(obj)
            s.coord  = obj.coord;
            s.connec = obj.connec;
            s.kFace  = obj.kFace -1;
            obj.mesh = Mesh.create(s);
        end
        
        function createGlobalConnec(obj)
            nodes = find(obj.nodesInBoxFaces);
            for inode = 1:size(obj.connec,2)
               obj.globalConnec(:,inode) = nodes(obj.connec(:,inode));
            end
        end
        
    end
    
end