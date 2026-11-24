classdef ProjectedGradientXiBounds < ProjectedGradient

    properties (Access = private)
        fixedUpperBound
        fixedLowerBound
        rhoIndices
        xiIndices
        mU
    end

    methods (Access = public)

        function obj = ProjectedGradientXiBounds(cParams)
            obj@ProjectedGradient(cParams);
            obj.fixedUpperBound = cParams.ub;
            obj.fixedLowerBound = cParams.lb;
            obj.mU              = cParams.mU;
            obj.createIndices(cParams.fieldSizes,cParams.rhoField,cParams.xiField);
        end

        function designVariable = update(obj,g,designVariable)
            y   = designVariable.fun.fValues - obj.tau*g;
            rho = obj.projectRho(y(obj.rhoIndices));
            [xiMin,xiMax] = XiLatticeDatasetGenerator.computeXiBounds(sqrt(1-rho),obj.mU);
            ub = obj.fixedUpperBound;
            lb = obj.fixedLowerBound;
            ub(obj.xiIndices) = xiMax;
            lb(obj.xiIndices) = xiMin;
            obj.updateBounds(ub,lb);
            designVariable = update@ProjectedGradient(obj,g,designVariable);
        end

    end

    methods (Access = private)

        function createIndices(obj,fieldSizes,rhoField,xiField)
            offsets        = [0,cumsum(fieldSizes(:)')];
            obj.rhoIndices = offsets(rhoField) + (1:fieldSizes(rhoField));
            obj.xiIndices  = offsets(xiField)  + (1:fieldSizes(xiField));
        end

        function rho = projectRho(obj,yRho)
            ub  = obj.fixedUpperBound(obj.rhoIndices);
            lb  = obj.fixedLowerBound(obj.rhoIndices);
            rho = min(ub,max(yRho,lb));
        end

    end

end