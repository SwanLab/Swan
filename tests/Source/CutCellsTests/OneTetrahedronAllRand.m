classdef OneTetrahedronAllRand <  VectorizedTriangulationTest
    
    properties (Access = private)
       testName = 'OneTetrahedronAllRand';        
    end
    
    methods (Access = protected)
        
        function createCoordAndConnec(obj)
            obj.coord = rand(4,3);              
            obj.connec = randperm(4);
        end
        
        function createLevelSet(obj)
            b = 10;
            a = 0;
            ls = rand(4,1);
            ls = a + (b-a)*ls;            
            difPos = randperm(4,1);
            isPositive = randperm(2,1);
            position = false(4,1);
            if isPositive == 1
                position(difPos) = true;
            else 
                position(setdiff(1:4,difPos)) = true;
            end            
            ls(position) = -ls(position);
            obj.levelSet = ls;
        end        
        
    end
    
end