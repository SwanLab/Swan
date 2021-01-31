classdef OneTetrahedronRandCoordOrderedConnecCase11 <  VectorizedTriangulationTest
    
    properties (Access = private)
       testName = 'OneTetrahedronRandCoordOrderedConnecCase11';        
    end
    
    methods (Access = protected)
        
        function init(obj)
            obj.coord = rand(4,3);            
            obj.connec = [1 2 3 4];
            obj.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];     
            obj.levelSet = [-7.8496;-9.7731;8.3404;-8.3622];            
        end
               
    end
    
end