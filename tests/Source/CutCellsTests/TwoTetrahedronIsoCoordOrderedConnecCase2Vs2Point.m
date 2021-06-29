classdef TwoTetrahedronIsoCoordOrderedConnecCase2Vs2Point <  VectorizedTriangulationTest
    
    properties (Access = private)
       testName = 'TwoTetrahedronIsoCoordOrderedConnecCase2Vs2Point';        
    end
    
    methods (Access = protected)
        
        function createCoordAndConnec(obj)
            obj.coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 1];            
            obj.connec = [1 2 3 4;
                          2 3 4 5];
            obj.boundaryConnec = [1 2 3;1 2 4;1 3 4;5 2 3;5 2 4;5 3 4];                        
        end
        
        function createLevelSet(obj)
             obj.levelSet = [7.8496;5.7731;-8.3404;-8.3622;10.3622];
        end        
        
    end
    
end