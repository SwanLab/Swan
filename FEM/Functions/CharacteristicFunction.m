classdef CharacteristicFunction < L2Function

    properties (Access = public)
        ndimf
    end

    properties (Access = private)
        levelSet
        backgroundMesh
        fValues
    end

    properties (Access = private)
        fieldFunction
    end

    methods (Access = public)

        function obj = CharacteristicFunction(cParams)
            obj.init(cParams);
            obj.computeFieldFunction();
        end

        function fxV = evaluate(obj,xV)
            s.fValues = obj.levelSet.value;
            s.mesh    = obj.backgroundMesh;
            ls        = P1Function(s);
            f         = ls.evaluate(xV);
            fxV       = obj.fieldFunction.evaluate(xV);
            nGaus     = size(f,2);
            for iGaus = 1:nGaus
                fG                = squeeze(f(1,iGaus,:));
                fxV(:,iGaus,fG>0) = 0;
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.ndimf          = 1;
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.levelSet.mesh;
            if isfield(cParams,'F')
                obj.fValues = cParams.F;
            else
                obj.fValues = ones(size(obj.levelSet.value));
            end
        end

        function computeFieldFunction(obj)
            s.fValues         = obj.fValues;
            s.mesh            = obj.backgroundMesh;
            obj.fieldFunction = P1Function(s);
        end

    end
end