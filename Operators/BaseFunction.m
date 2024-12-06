classdef BaseFunction < handle & matlab.mixin.Copyable

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (GetAccess = public, SetAccess = protected)
        mesh
        ndimf
    end

    methods (Access = public)

        function fun = project(obj,target)
            s.mesh          = obj.mesh;
            s.projectorType = target;
            proj = Projector.create(s);
            fun = proj.project(obj);
        end       

        function plot(obj)
            p1D = project(obj,'P1D');
            p1D.plot();
        end    

        function plotVector(obj,varargin)
            if size(varargin, 1) == 1, n = varargin{1}; else, n = 2; end
            p1D = project(obj,'P1D');            
            plotVector(p1D,n);
        end             

        function r = ctranspose(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) pagetranspose(aOp(xV));
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = plus(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) + bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = minus(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) - bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = times(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            ndimfA = BaseFunction.computeFieldDimension(a);
            ndimfB = BaseFunction.computeFieldDimension(b);
            s.operation = @(xV) aOp(xV).*bOp(xV);
            s.ndimf = max(ndimfA,ndimfB);
            if isa(a,'BaseFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = mtimes(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) pagemtimes(aOp(xV),bOp(xV));
            if isa(a,'BaseFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = rdivide(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV)./bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = uminus(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) -aOp(xV);
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = power(a,b)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) aOp(xV).^b;
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = norm(a,b)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) pagenorm(aOp(xV),b);
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = log(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) log(aOp(xV));
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = exp(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) exp(aOp(xV));
            s.mesh = a.mesh;
            r = DomainFunction(s);
        end

        function r = trace(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) trace(aOp(xV));
            s.ndimf = a.ndimf;
            s.mesh  = a.mesh;
            r = DomainFunction(s);
        end

        function f = mrdivide(f1,f2)
            s.operation = @(xV) f1.evaluate(xV)./f2.evaluate(xV);
            s.ndimf = max(f1.ndimf,f2.ndimf);
            s.mesh = f1.mesh;
            f = DomainFunction(s);
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

    methods (Access = public, Abstract)
        evaluate(obj,xV)
    end

end