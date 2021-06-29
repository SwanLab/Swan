classdef FindingBug <  VectorizedTriangulationTest
    
    properties (Access = private)
       testName = 'FindingBug';        
    end
    
    methods (Access = protected)
        
        function createCoordAndConnec(obj)
            obj.coord = [0      0        0;
                         1      0       0;                         
                         0      1        0;                         
                         0.1   0.1     0.7];
            obj.connec = [2     1     3     4];
            obj.boundaryConnec = [2 4 1;2 4 3;2 1 3;4 1 3];      
        end
        
        function createLevelSet(obj)
             obj.levelSet = [-3.9507;15.2486;-10.0641;5.1921];
        end        
        
    end
    
end