classdef HessianComputer < handle

    properties (Access = private)
        initHessian
        designVariable
        cost
    end

    properties (Access = private)
        deltaX
        deltaCostGrad
    end

    methods (Access = public)
        function obj = HessianComputer(cParams)
            obj.init(cParams)
        end

        function newHessian = compute(obj)
            obj.computeDesignVariableDiff();
            obj.computeCostGradientDiff();
            newHessian = obj.computeNewHessian();
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.initHessian    = cParams.initHessian;
            obj.designVariable = cParams.designVariable;
            obj.cost           = cParams.cost;
        end

        function computeDesignVariableDiff(obj)
            x          = obj.designVariable.value;
            xOld       = obj.designVariable.valueOld;
            obj.deltaX = x-xOld;
        end

        function computeCostGradientDiff(obj)
            DJ                = obj.cost.gradient;
            DJOld             = obj.cost.gradientOld;
            obj.deltaCostGrad = DJ-DJOld;
        end

        function newHessian = computeNewHessian(obj)
            dimX = length(obj.designVariable.value);
            H    = obj.initHessian;
            dX   = obj.deltaX;
            dDJ  = obj.deltaCostGrad;
            if dDJ == zeros(size(dDJ))
                newHessian = zeros(dimX,dimX);
            else
                costFrac   = (dDJ*dDJ')/(dDJ'*dX);
                xFrac      = (H*dX*(H*dX)')/(dX'*H*dX);
                newHessian = H + costFrac - xFrac;
            end
        end
    end
end