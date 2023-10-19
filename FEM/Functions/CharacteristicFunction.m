classdef CharacteristicFunction < L2Function

    properties (Access = public)
        ndimf
    end

    properties (Access = private)
        levelSet
        backgroundMesh
    end

    methods (Access = public)

        function obj = CharacteristicFunction(cParams)
            obj.init(cParams);
        end

        function fxV = evaluate(obj,xV)
            s.fValues = obj.levelSet.value;
            s.mesh    = obj.backgroundMesh;
            ls        = P1Function(s);
            f         = ls.evaluate(xV);
            nGaus     = size(f,2);
            nElem     = size(f,3);
            fxV       = zeros(1,nGaus,nElem);
            for iGaus = 1:nGaus
                fG                 = squeeze(f(1,iGaus,:));
                fxV(1,iGaus,fG<=0) = 1;
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ndimf          = 1;
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.levelSet.mesh;
        end

    end
end