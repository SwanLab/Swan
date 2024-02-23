classdef MathEntity < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = protected)
        
        function obj = MathEntity(cParams)
        end
              
        
        function  power(b)
            obj = @(x) x.^b;
        end
    end
    
end