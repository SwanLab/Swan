classdef OneTetrahedronAllRand2Vs2 <  VectorizedTriangulationTest
    
    properties (Access = private)
       testName = 'OneTetrahedronAllRand2Vs2';        
    end
    
    methods (Access = protected)
        
        function createCoordAndConnec(obj)
            obj.coord = rand(4,3);              
            obj.connec = randperm(4);
            obj.boundaryConnec = [1 2 3;1 2 4;1 3 4;2 3 4];                        
        end
        
        function createLevelSet(obj)
            b = 10;
            a = 0;
            ls = rand(4,1);
            ls = a + (b-a)*ls;            
            position = randperm(4,2);      
            ls(position) = -ls(position);
            obj.levelSet = ls;
        end        
        
    end
    
end