classdef BaseFunction < handle & matlab.mixin.Copyable

    properties (Access = public)

    end

    properties (GetAccess = protected, SetAccess = protected)
       fxVOld
       xVOld
    end

    properties (GetAccess = public, SetAccess = protected)
        mesh
        ndimf
        ndimfTotal
    end

    methods (Access = public)

        function fxV = evaluate(obj, xV)
            if ~isequal(xV,obj.xVOld) || isempty(obj.fxVOld)
                fxV = obj.evaluateNew(xV);
                obj.fxVOld  = fxV;
                obj.xVOld = xV;
            else
                fxV = obj.fxVOld;
            end
        end

        function fun = project(obj,target)
            switch class(obj) % Parche 1: ndimF
                case {'UnfittedFunction','UnfittedBoundaryFunction'}
                    nTensor = 1;
                otherwise
                    nTensor = length(obj.ndimf);
            end

            if nTensor>=4 % Parche 2
                s.projectorType = target;
                proj = ProjectorToLagrangianTensor(s);
                fun = proj.project(obj);
            else
                s.mesh          = obj.mesh;
                s.projectorType = target;
                proj = Projector.create(s);
                fun = proj.project(obj);
            end
        end

        function plot(obj)
            p1D = project(obj,'P1D');
            p1D.plot();
        end    

        function print(obj,varargin)
            p1D = project(obj,'P1D');
            p1D.print(varargin{:});
        end

        function plotVector(obj,varargin)
            if size(varargin, 1) == 1, n = varargin{1}; else, n = 2; end
            p1D = project(obj,'P1D');            
            plotVector(p1D,n);
        end             

        function r = ctranspose(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) pagetranspose(aOp(xV));
            s.mesh  = a.mesh;
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function r = plus(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) + bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh  = a.mesh;
                s.ndimf = a.ndimf; 
            else
                s.mesh  = b.mesh;
                s.ndimf = b.ndimf; 
            end
            r = DomainFunction(s);
        end

        function r = minus(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) - bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh  = a.mesh;
                s.ndimf = a.ndimf; 
            else
                s.mesh  = b.mesh;
                s.ndimf = b.ndimf; 
            end            
            r = DomainFunction(s);
        end

        function r = times(a,b)
            ndimfA = BaseFunction.computeFieldDimension(a);
            ndimfB = BaseFunction.computeFieldDimension(b);
            nTensorA = length(ndimfA);
            nTensorB = length(ndimfB);
            a = Expand(a,max(nTensorB,2)); b = Expand(b,max(nTensorA,2));
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) squeezeParticular(aOp(xV).*bOp(xV),2);
            s.ndimf = max(ndimfA,ndimfB);
            if isa(a,'BaseFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        function r = mtimes(a,b)
            ndimfA = BaseFunction.computeFieldDimension(a);
            ndimfB = BaseFunction.computeFieldDimension(b);
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) pagemtimes(aOp(xV),bOp(xV));
            s.ndimf = [ndimfA(1) ndimfB(1)];
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
            ndimfA = BaseFunction.computeFieldDimension(a);
            ndimfB = BaseFunction.computeFieldDimension(b);
            s.operation = @(xV) aOp(xV)./bOp(xV);
            s.ndimf = max(ndimfA,ndimfB);
            if isa(a,'BaseFunction')
                s.mesh = a.mesh;
            else
                s.mesh = b.mesh;
            end
            r = DomainFunction(s);
        end

        % function r = mrdivide(a,b)
        %     r = obj.rdivide(a,b);
        % end

        function r = uminus(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) -aOp(xV);
            s.mesh = a.mesh;
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function r = power(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV).^bOp(xV);
            s.mesh  = a.mesh;
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function r = sqrt(a)
            r = power(a,0.5);
        end

        function r = gt(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) > bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh  = a.mesh;
                s.ndimf = a.ndimf; 
            else
                s.mesh  = b.mesh;
                s.ndimf = b.ndimf; 
            end  
            r = DomainFunction(s);
        end

        function r = ge(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) >= bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh  = a.mesh;
                s.ndimf = a.ndimf; 
            else
                s.mesh  = b.mesh;
                s.ndimf = b.ndimf; 
            end  
            r = DomainFunction(s);
        end
        function r = lt(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) < bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh  = a.mesh;
                s.ndimf = a.ndimf; 
            else
                s.mesh  = b.mesh;
                s.ndimf = b.ndimf; 
            end  
            r = DomainFunction(s);
        end

        function r = le(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) aOp(xV) <= bOp(xV);
            if isa(a,'BaseFunction')
                s.mesh  = a.mesh;
                s.ndimf = a.ndimf; 
            else
                s.mesh  = b.mesh;
                s.ndimf = b.ndimf; 
            end  
            r = DomainFunction(s);
        end

        function r = not(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) ~aOp(xV);
            s.mesh  = a.mesh;
            s.ndimf = a.ndimf;
            r = DomainFunction(s);
        end

        function r = norm(varargin)
            a = varargin{1};
            if nargin == 1
                b = 2;
            elseif nargin == 2
                b = varargin{2};
            end
            a = Expand(a); % To perform pagenorm of a vector (ndim -> ndim x 1)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) squeezeParticular(pagenorm(aOp(xV),b),2);
            s.mesh = a.mesh;
            s.ndimf = 1;            
            r = DomainFunction(s);
        end

        function r = log(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) log(aOp(xV));
            s.mesh  = a.mesh;
            s.ndimf = a.ndimf;            
            r = DomainFunction(s);
        end

        function r = exp(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) exp(aOp(xV));
            s.mesh  = a.mesh;
            s.ndimf = a.ndimf;            
            r = DomainFunction(s);
        end

        function r = abs(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) abs(aOp(xV));
            s.mesh  = a.mesh;
            s.ndimf = a.ndimf;            
            r = DomainFunction(s);
        end        

        function r = trace(a)
            aOp = BaseFunction.computeOperation(a);
            s.operation = @(xV) trace(aOp(xV));
            s.ndimf = 1;
            s.mesh  = a.mesh;
            r = DomainFunction(s);
        end

        function r = min(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) min(aOp(xV),bOp(xV));
            if isa(a,'BaseFunction')
                s.mesh  = a.mesh;
                s.ndimf = a.ndimf; 
            else
                s.mesh  = b.mesh;
                s.ndimf = b.ndimf; 
            end  
            r = DomainFunction(s);
        end

        function r = max(a,b)
            aOp = BaseFunction.computeOperation(a);
            bOp = BaseFunction.computeOperation(b);
            s.operation = @(xV) max(aOp(xV),bOp(xV));
            if isa(a,'BaseFunction')
                s.mesh  = a.mesh;
                s.ndimf = a.ndimf; 
            else
                s.mesh  = b.mesh;
                s.ndimf = b.ndimf; 
            end  
            r = DomainFunction(s);
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
            if isnumeric(a)
                if isscalar(a)
                    ndimf = 1;
                else
                    ndimf = size(a);
                end
            else
                ndimf = a.ndimf;
            end
        end

    end

    methods (Access = protected, Abstract)
        evaluateNew(obj,xV)
    end

end