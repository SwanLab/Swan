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

        function r = times(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV).*bOp(xV);
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

    end

    methods (Access = public, Static)
        
        function obj = create(handle,l2F)
            s.handleFunction = handle;
            s.l2function     = l2F;
            obj = CompositionFunction(s);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.handleFunction = cParams.handleFunction;
            obj.l2function     = cParams.l2function;
            obj.ndimf = size(obj.handleFunction,2);
        end
        
    end
    
end