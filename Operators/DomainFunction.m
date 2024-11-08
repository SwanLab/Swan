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

        function plot(obj,m)
            fD = obj.project('P1D',m);
            fD.plot();
        end

        function r = ctranspose(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) nOrderTranspose(aOp(xV));
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
            if not(isfloat(a))
                aOp = DomainFunction.computeOperation(a);
            else
                aOp = @(xV) a;
            end
            if not(isfloat(b))
                bOp = DomainFunction.computeOperation(b);
            else
                bOp = @(xV) b;
            end
            s.operation = @(xV) aOp(xV).*bOp(xV);
            r = DomainFunction(s);
        end

        function r = mtimes(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) pagemtimes(aOp(xV),bOp(xV));
            r = DomainFunction(s);
        end

        function r = rdivide(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV)./bOp(xV);
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

        function r = exp(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) exp(aOp(xV));
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
            if isfield(cParams,'ndimf')
                obj.ndimf = cParams.ndimf;
            else
                obj.ndimf = 1;
            end
        end

    end

    methods (Static, Access = public)

        function op = computeOperation(a)
            if isprop(a,'operation')
                op = a.operation;
            elseif isnumeric(a)
                op = @(xV) a;
            else
                op = @(xV) a.evaluate(xV);
            end
        end

    end

end