classdef DomainFunction < handle
    
    properties (Access = public)
        operation
    end
    
    properties (Access = private)
        
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
        
        function r = ctranspose(a)
            s.operation = @(xV) pagetranspose(a.operation(xV));
            r = DomainFunction(s);
        end
        
        function r = plus(a,b)
            s.operation = @(xV) a.operation(xV) + b.operation(xV);
            r = DomainFunction(s);
        end
        
        function r = minus(a,b)
            s.operation = @(xV) a.operation(xV) - b.operation(xV);
            r = DomainFunction(s);
        end
        
        function r = uminus(a)
            s.operation = @(xV) -a.operation(xV);
            r = DomainFunction(s);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.operation = cParams.operation;
        end
        
    end
    
end