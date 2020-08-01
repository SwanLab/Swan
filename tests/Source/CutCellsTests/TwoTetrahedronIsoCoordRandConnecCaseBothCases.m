classdef TwoTetrahedronIsoCoordRandConnecCaseBothCases <  VectorizedTriangulationTest
    
    properties (Access = private)
       testName = 'TwoTetrahedronIsoCoordRandConnecCaseBothCases';        
    end
    
    methods (Access = protected)
        
        function createCoordAndConnec(obj)
            obj.coord = rand(50,3); 
            t = delaunayTriangulation(obj.coord);
            obj.connec = t.ConnectivityList;
            obj.boundaryConnec = boundary(obj.coord);                        
        end
        
        function createLevelSet(obj)
            b = 10;
            a = 0;
            ls = rand(size(obj.coord,1),1);
            ls = a + (b-a)*ls;            
            number = randperm(size(obj.coord,1),1);
            position = randperm(size(obj.coord,1),number);      
            ls(position) = -ls(position);
            obj.levelSet = ls;
        end        
        
    end
    
end