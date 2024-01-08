classdef CompositionFunction < L2Function
    
    properties (Access = public)
        ndimf
    end
    
    properties (Access = private)
        handleFunction
        l2function
    end
    
    methods (Access = public)
        
        function obj = CompositionFunction(cParams)
            obj.init(cParams)            
        end

        function fV = evaluate(obj,xV)
            l2fV = obj.l2function.evaluate(xV);
            fV   = obj.handleFunction(l2fV);
        end

        function plot(obj)
            f = obj.project('P1D');
            f.plot();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.handleFunction = cParams.handleFunction;
            obj.l2function     = cParams.l2function;
            obj.ndimf = size(obj.handleFunction,2);
        end
        
    end
    
end