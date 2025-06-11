classdef AnalyticalFunction2 < BaseFunction

    properties (Access = private)
        domainFunction
    end

    methods (Access = public)

        function obj = AnalyticalFunction2(cParams)
            obj.init(cParams)
        end

    end

    methods (Access = public, Static)

        function obj = create(fHandle,mesh)
            s.fHandle = fHandle;
            s.mesh    = mesh;
            obj = AnalyticalFunction2(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            fHandle   = cParams.fHandle;
            obj.mesh  = cParams.mesh;
            obj.computeNdimf(fHandle);
            obj.createDomainFunction(fHandle);
        end

        function computeNdimf(obj,fHandle)
            ndim = obj.mesh.ndim;
            obj.ndimf = size(fHandle(zeros(ndim,1)),1);
        end

        function createDomainFunction(obj,fHandle)
            s.operation = @(xV) obj.compute(fHandle,xV);
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
