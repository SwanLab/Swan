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
            LamDirection = [0 1 0];
            obj.LaminateDirection = LamDirection/norm(LamDirection);
        end
        
        function loadFiberDirection(obj)
            Direction = [1 0 0];   
            obj.FiberDirection = Direction/norm(Direction);
        end
        
    end
   
end