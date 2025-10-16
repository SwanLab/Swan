classdef AnalyticalFunction < BaseFunction

    properties (Access = private)
        fHandle
        domainFunction
    end

    methods (Access = public)

        function obj = AnalyticalFunction(cParams)
            obj.init(cParams)
        end

        function f = createNew(obj,mesh)
            f = AnalyticalFunction.create(obj.fHandle,mesh);
        end

    end

    methods (Access = public, Static)

        function obj = create(fHandle,mesh)
            s.fHandle = fHandle;
            s.mesh    = mesh;
            obj = AnalyticalFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fHandle = cParams.fHandle;
            obj.mesh    = cParams.mesh;
            obj.computeNdimf();
            obj.createDomainFunction();
        end

        function computeNdimf(obj)
            ndim = obj.mesh.ndim;
            obj.ndimf = size(obj.fHandle(zeros(ndim,1)));
            obj.ndimfTotal = prod(obj.ndimf);
        end

        function createDomainFunction(obj)
            s.operation = @(xV) obj.compute(obj.fHandle,xV);
            s.mesh      = obj.mesh;
            s.ndimf     = obj.ndimf;
            f = DomainFunction(s);
            obj.domainFunction = f;
        end

        function fxV = compute(obj,fHandle, xV)
            x   = obj.mesh.computeXgauss(xV);
            fxV = fHandle(x);
        end

    end

    methods (Access = protected)

        function fxV = evaluateNew(obj, xV)
            fxV = obj.domainFunction.evaluate(xV);
        end

    end

end
