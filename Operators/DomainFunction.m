classdef DomainFunction < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        operation
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = DomainFunction(cParams)
            obj.init(cParams)
            
        end
        
        function r = evaluate(obj,xV)
            r = obj.operation(xV);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.operation = cParams.operation;
        end
        
    end
    
end