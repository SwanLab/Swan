classdef HomogenizedMicrostructureInterpolator < handle

    properties (Access = private)
        density
        fileName
        Chomog
        mesh
    end

    methods (Access = public)

        function obj = HomogenizedMicrostructureInterpolator(cParams)
            obj.init(cParams)
        end

        function C = obtainTensor(obj)
            fun = obj.Chomog.fun;
            s.operation = @(xV) obj.evaluate(fun,xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj)
            fun = obj.Chomog.dfun;
            s.operation = @(xV) obj.evaluate(fun,xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            dC{1} =  DomainFunction(s);
        end

        function d2C = obtainTensorSecondDerivative(obj)
            fun = obj.Chomog.ddfun;
            s.operation = @(xV) obj.evaluate(fun,xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            d2C{1} =  DomainFunction(s);
        end

        function setDesignVariable(obj,x)
            obj.density = x;            
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh   = cParams.mesh;
            obj.Chomog = cParams.Chomog;
        end

        function C = evaluate(obj,fun,xV)
            nStre = size(fun,1);
            nGaus = size(xV,2);
            nElem = obj.mesh.nelem;
            C = zeros(nStre,nStre,nStre,nStre,nGaus,nElem);
            phiV = obj.density{1}.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    for k=1:nStre
                        for l=1:nStre
                            C(i,j,k,l,:,:) = fun{i,j,k,l}(phiV);
                        end
                    end
                end
            end
        end

    end

end