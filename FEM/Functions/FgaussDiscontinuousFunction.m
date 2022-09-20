classdef FgaussDiscontinuousFunction < handle
    
    properties (Access = public)
        fValues
        
    end
    
    properties (Access = private)
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = FgaussDiscontinuousFunction(cParams)
            obj.init(cParams)
            
        end
        
        function fxV = evaluate(obj, xV)
            fxV = obj.fValues;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fValues = cParams.fValues;
            obj.xGauss  = cParams.xGauss;
        end
        
    end
    
end