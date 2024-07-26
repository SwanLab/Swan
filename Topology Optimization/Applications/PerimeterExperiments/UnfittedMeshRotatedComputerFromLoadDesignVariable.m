classdef UnfittedMeshRotatedComputerFromLoadDesignVariable < handle
    
    properties (Access = private)
        filePath
    end
    
    methods (Access = public)
        
        function obj = UnfittedMeshRotatedComputerFromLoadDesignVariable(cParams)
            obj.init(cParams)
        end
        
        function u = computeUnfittedMesh(obj)
            d        = load(obj.filePath);
            bMesh    = obj.turnMesh(d.mesh);
            levelSet = d.x;
            s.backgroundMesh = bMesh;
            s.boundaryMesh = obj.computeBoundaryMesh(bMesh);
            u = UnfittedMesh(s);
            u.compute(levelSet)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filePath = cParams.filePath;
        end
  
    end
    
    methods (Access = private, Static)
        
        function bMesh = computeBoundaryMesh(mesh)
            sB.backgroundMesh = mesh;
            sB.dimension = 1:3;
            sB.type = 'FromReactangularBox';
            bM = BoundaryMeshCreator.create(sB);
            bMesh = bM.create();
        end
   
        function m = turnMesh(mesh)
            s.coord = mesh.coord(:,[1 3 2]);
            s.connec = mesh.connec;
            m = Mesh.create(s);
        end
        
    end
    
end