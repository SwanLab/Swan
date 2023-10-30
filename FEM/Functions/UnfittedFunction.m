classdef UnfittedFunction < handle

    properties (Access = public)
        ndimf
        levelSet
        fun
    end

    properties (Access = private)
        backgroundMesh
    end

    methods (Access = public)

        function obj = UnfittedFunction(cParams)
            obj.init(cParams);
        end

        function fxV = evaluate(obj,xV)
            s.fValues = obj.levelSet.value;
            s.mesh    = obj.backgroundMesh;
            ls        = P1Function(s);
            f         = ls.evaluate(xV);
            fxV       = obj.fun.evaluate(xV);
            nGaus     = size(f,2);
            for iGaus = 1:nGaus
                fG                = squeeze(f(1,iGaus,:));
                fxV(:,iGaus,fG>0) = 0;
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.levelSet       = cParams.levelSet;
            obj.backgroundMesh = cParams.levelSet.mesh;
            obj.fun            = cParams.fun;
            obj.ndimf          = cParams.fun.ndimf;
        end

    end
end