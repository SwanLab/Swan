classdef AcceleratedMMA < handle

    properties (Access = public)
        xFinal
        incJ
        J
    end

    properties (Access = private)
        cost
        TOL
        maxIter
        momentumParameter
        nIter
        oldCost
        designVariable
        MMA
    end

    methods (Access = public)

        function obj = AcceleratedMMA(cParams)
            obj.init(cParams);
        end

        function solve(obj)
            d = obj.designVariable;
            obj.restartDesignVariable();
            obj.cost.computeFunctionAndGradient(d);
            x    = obj.designVariable.fun.fValues;
            xOld = x;
            obj.oldCost = 0;
            obj.nIter   = 1;
            obj.storeIterationData();
            while ~obj.hasConverged()
                beta = obj.computeBeta();
                y    = obj.addMomentumTerm(xOld,x,beta);
                xNew = obj.MMA.computeIteration(y);
                xOld = x;
                x    = xNew;
                d.update(x);
                obj.cost.computeFunctionAndGradient(d);
                obj.storeIterationData();
                obj.oldCost = obj.cost.value;
                obj.nIter   = obj.nIter + 1;
            end
            obj.xFinal = x;
        end

    end

    methods (Access = private, Static)

        function y = addMomentumTerm(xOld,xNew,beta)
            y = xNew + beta*(xNew - xOld);
        end

        function x = computeGradientStep(x,dJ,t)
            x = x - t*dJ;
        end

        function x = computeProjection(x)
            x = max(min(x,1),0);
        end

        function [xOld,xNew] = updateXnewXold(xNew,xNewNew)
            xOld = xNew;
            xNew = xNewNew;
        end

    end

    methods (Access = private)

        function beta = computeBeta(obj)
            s.iter = obj.nIter;
            beta   = obj.momentumParameter.computeValue(s);
        end

        function init(obj,cParams)
            obj.cost              = cParams.cost;
            obj.designVariable    = cParams.designVariable;            
            obj.TOL               = cParams.TOL;
            obj.maxIter           = cParams.maxIter;
            obj.momentumParameter = MomentumParameter.create(cParams.momentumParameter);
            obj.MMA               = Nesterov_MMA(cParams);
        end

        function storeIterationData(obj)
            obj.incJ(obj.nIter) = abs(obj.cost.value - obj.oldCost);
            obj.J(obj.nIter)    = obj.cost.value;
        end

        function restartDesignVariable(obj)
            obj.designVariable.update(ones(numel(obj.designVariable.fun.fValues),1));
        end

        function itHas = hasConverged(obj)
            if obj.nIter == 1
                itHas = false;
            else
                itHas = obj.nIter >= obj.maxIter || obj.incJ(end) < obj.TOL;
            end
        end

    end

end

