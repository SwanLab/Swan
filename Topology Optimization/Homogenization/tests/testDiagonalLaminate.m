classdef testDiagonalLaminate < TestSequentialLaminateTestedWithNumerics

    
    properties 
    end
    
    methods (Access = public)
        
        function obj = testDiagonalLaminate()            
            obj.compute()            
        end
        
    end
    
    methods (Access = protected)
        
        function loadLaminateDirection(obj)
            Direction = [-1 1 0];
            obj.LaminateDirection = Direction/norm(Direction);
        end
        
        function loadFiberDirection(obj)
            Direction = [1 1 0];   
            obj.FiberDirection = Direction/norm(Direction);
        end
        
    end
   
end