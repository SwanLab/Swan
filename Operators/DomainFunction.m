classdef DomainFunction < handle

    properties (Access = public)
        ndimf
    end

    properties (GetAccess = public, SetAccess = protected)
        operation
        mesh
    end
    
    methods (Access = public)
        
        function obj = DomainFunction(cParams)
            obj.init(cParams)
        end
        
        function r = evaluate(obj,xV)
            r = obj.operation(xV);
        end

        function plot(obj)
            fD = obj.project('P1D',obj.mesh);
            fD.plot();
        end

        function fun = project(obj,target,mesh)
            s.mesh          = mesh;
            s.projectorType = target;
            proj = Projector.create(s);
            fun = proj.project(obj);
        end

        function r = ctranspose(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) pagetranspose(aOp(xV));
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end
        
        function r = plus(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) + bOp(xV);
            if isa(a,'DomainFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = minus(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) - bOp(xV);
            if isa(a,'DomainFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = times(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            ndimfA = DomainFunction.computeFieldDimension(a);
            ndimfB = DomainFunction.computeFieldDimension(b);
            s.operation = @(xV) aOp(xV).*bOp(xV);
            s.ndimf = max(ndimfA,ndimfB);
            if isa(a,'DomainFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = mtimes(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) pagemtimes(aOp(xV),bOp(xV));
            if isa(a,'DomainFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = rdivide(a,b)
            aOp = DomainFunction.computeOperation(a);
            bOp = DomainFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV)./bOp(xV);
            if isa(a,'DomainFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end
        
        function r = uminus(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) -aOp(xV);
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = power(a,b)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) aOp(xV).^b;
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = norm(a,b)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) pagenorm(aOp(xV),b);
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = log(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) log(aOp(xV));
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = exp(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) exp(aOp(xV));
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end        

        function r = trace(a)
            aOp = DomainFunction.computeOperation(a);
            s.operation = @(xV) trace(aOp(xV));
            s.ndimf = a.ndimf;
            s.mesh = a.mesh;
            r = DomainFunction(s);
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

            if isfield(cParams,'mesh')
                obj.mesh = cParams.mesh;
            else
                obj.mesh = [];
            end
        end

    end

    methods (Static, Access = public)

        function op = computeOperation(a)
            if isa(a,'DomainFunction')
                op = a.operation;
            elseif isnumeric(a)
                op = @(xV) a;
            else
                op = @(xV) a.evaluate(xV);
            end
        end

        function ndimf = computeFieldDimension(a)
            if isa(a,'DomainFunction')
                ndimf = a.ndimf;
            elseif isnumeric(a)
                ndimf = size(a,1);
            else
                ndimf = a.ndimf;
            end
        end        



    end

end