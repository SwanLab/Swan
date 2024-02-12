classdef ShFunc_ExternalWork < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ShFunc_ExternalWork(cParams)
            obj.init(cParams)
            
        end
        
        function F = computeFunction(obj)
        end
        
        function J = computeGradient(obj)
        end
        
        function H = computeHessian(obj)
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            
        end
        
    end
    
end