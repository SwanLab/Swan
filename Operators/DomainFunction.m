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
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) pagetranspose(aOp(xV));
            r = DomainFunction(s);
        end
        
        function r = plus(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) + bOp(xV);
            r = DomainFunction(s);
        end

        function r = minus(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) - bOp(xV);
            r = DomainFunction(s);
        end

        function r = times(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV).*bOp(xV);
            r = DomainFunction(s);
        end
        
        function r = uminus(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) -aOp(xV);
            r = DomainFunction(s);
        end

        function r = power(a,b)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) aOp(xV).^b;
            r = DomainFunction(s);
        end

        function r = norm(a,b)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) pagenorm(aOp(xV),b);
            r = DomainFunction(s);
        end

        function r = log(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) log(aOp(xV));
            r = DomainFunction(s);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.operation = cParams.operation;
        end
        
    end

    methods (Static, Access = public)

        function op = computeOperation(a)
            if isprop(a,'operation')
                op = a.operation;
            else
                op = @(xV) evaluate(xV);
            end
        end

    end

end