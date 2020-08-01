classdef VectorizedTriangulationTest < testShowingError
    
    properties (Access = protected)
        tol = 1e-14;
    end

    properties (Access = private)
        backgroundMesh                
    end
    
    properties (Access = protected)
        coord          
        connec
        boundaryConnec        
        levelSet        
    end
    
    methods (Access = public)
        
        function obj = VectorizedTriangulationTest()
            obj.createCoordAndConnec();
            obj.createLevelSet();
            obj.createBackgroundMesh();             
        end
        
    end
    
    methods (Access = private)

        function createBackgroundMesh(obj)
            s.connec = obj.connec;
            s.coord  = obj.coord;
            m = Mesh(s);
            obj.backgroundMesh = m;
        end                
        
    end
    
    methods (Access = protected)
       
        function computeError(obj)
            s.backgroundMesh = obj.backgroundMesh;
            s.levelSet       = obj.levelSet;
            s.boundaryConnec = obj.boundaryConnec;
            c = ComputingCutMeshVectorized(s);
            e = c.compute();
            obj.error = e;
        end                 
        
    end
    
    methods (Access = protected, Abstract)
        createCoordAndConnec(obj)
        createLevelSet(obj)
    end
    
end