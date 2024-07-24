classdef DomainFunction < handle
    
    properties (Access = public)
        operation
        ndimf
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
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end
        
        function r = plus(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) + bOp(xV);
            s.ndimf = b.ndimf;
            r = DomainFunction(s);
        end

        function r = minus(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) - bOp(xV);
            s.ndimf = b.ndimf;
            r = DomainFunction(s);
        end

        function r = times(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV).*bOp(xV);
            s.ndimf = b.ndimf;
            r = DomainFunction(s);
        end
        
        function r = uminus(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) -aOp(xV);
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function r = power(a,b)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) aOp(xV).^b;
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function r = norm(a,b)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) squeezeParticular(pagenorm(aOp(xV),b),1);
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function r = log(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) log(aOp(xV));
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function r = trace(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) trace(aOp(xV));
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function fun = project(obj,target,mesh)
            s.mesh          = mesh;
            s.projectorType = target;
            proj = Projector.create(s);
            fun = proj.project(obj);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.operation = cParams.operation;
            obj.ndimf = cParams.ndimf;
        end
        
    end

    methods (Static, Access = public)

        function op = computeOperation(a)
            if isprop(a,'operation')
                op = a.operation;
            elseif class(a) == "double"
                op = @(xV) a;
            else
                op = @(xV) a.evaluate(xV);
            end
        end

    end

end