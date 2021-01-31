classdef OneTetrahedronIsoCoordOrderedConnecCase7 <  VectorizedTriangulationTest
    
    properties (Access = private)
       testName = 'OneTetrahedronIsoCoordOrderedConnecCase7';        
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.coord    = [0 0 0; 1 0 0; 0 1 0; 0 0 1];            
            obj.connec   = [1 2 3 4];
            obj.levelSet = [-7.8496;-9.7731;-8.3404;8.3622];            
            obj.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4]; 
            obj.connecBcutMesh = [5 6 7];            
        end
    
    end
    
end