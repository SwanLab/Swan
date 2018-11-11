classdef testHorizontalLaminate < TestSequentialLaminateTestedWithNumerics

    
    properties 
    end
    
    methods (Access = public)
        
        function obj = testHorizontalLaminate()            
            obj.compute()            
        end
        
    end
    
    methods (Access = protected)
        
        function loadLaminateDirection(obj)
            obj.LaminateDirection = [0 1 0];            
        end
        
        function loadFiberDirection(obj)
            obj.FiberDirection = [1 0 0];            
        end
        
    end
   
end