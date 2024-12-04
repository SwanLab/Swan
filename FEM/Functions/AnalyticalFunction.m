classdef AnalyticalFunction < BaseFunction

    properties (Access = private)
        domainFunction
    end

    methods (Access = public)

        function obj = AnalyticalFunction(cParams)
            obj.init(cParams)
        end

        function fxV = evaluate(obj, xV)
            fxV = obj.domainFunction.evaluate(xV);
        end

    end

    methods (Access = public, Static)

        function obj = create(fHandle,ndimf,mesh)
            s.fHandle = fHandle;
            s.ndimf   = ndimf;
            s.mesh    = mesh;
            obj = AnalyticalFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            fHandle   = cParams.fHandle;
            obj.ndimf = cParams.ndimf;
            obj.mesh  = cParams.mesh;
            obj.createDomainFunction(fHandle);
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


end
