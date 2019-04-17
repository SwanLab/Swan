classdef testComputingTopOptWithVademecumData < testShowingError ...
   
    
    properties (Access = protected)
        tol = 1e-6;
        testName = 'testComputingTopOptWithVademecumData';        
    end
    
    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = testComputingTopOptWithVademecumData()
           
        end
    end
    
    methods (Access = protected)
        
        function computeError(obj)
            obj.error = 0;     
        end        
        
        
    end
    
    methods (Access = private)
        
   
        
    end
    
end