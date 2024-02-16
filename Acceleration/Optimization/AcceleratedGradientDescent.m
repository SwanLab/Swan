classdef AcceleratedGradientDescent < handle

    properties (Access = public)
        xFinal
        incJ
        J
    end

    properties (Access = private)
        cost
        tau
        TOL
        maxIter
        momentumParameter
        gradientComputer
        nIter
        oldCost
        designVariable
    end

    methods (Access = public)

        function obj = AcceleratedGradientDescent(cParams)
            obj.init(cParams);
        end

        function solve(obj)
            d = obj.designVariable;
            obj.restartDesignVariable();
            x    = obj.designVariable.fun.fValues;
            xOld = x;
            obj.cost.computeFunctionAndGradient(d);
            obj.oldCost = 0;
            obj.nIter   = 1;
            obj.storeIterationData();
            t           = obj.tau;
            while ~obj.hasConverged()
                beta = obj.computeBeta();
                y    = obj.addMomentumTerm(xOld,x,beta);
                obj.gradientComputer(x,y);
                DJ   = obj.cost.gradient;
                xNew = obj.computeGradientStep(y,DJ,t);
                xNew = obj.computeProjection(xNew);
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
            obj.tau               = cParams.tau;
            obj.momentumParameter = MomentumParameter.create(cParams.momentumParameter);
            obj.gradientComputer  = obj.createGradient(cParams.gDescentType);
        end

        function storeIterationData(obj)
            obj.incJ(obj.nIter) = abs(obj.cost.value - obj.oldCost);
            obj.J(obj.nIter)    = obj.cost.value;
        end

        function f = createGradient(obj,type)
            switch type
                case 'Polyak'
                    f = @(x,y) obj.computePolyakGradient(x,y);
                case 'Nesterov'
                    f = @(x,y) obj.computeNesterovGradient(x,y);
                otherwise
                    error('Case not implemented.')
            end
        end

        function computePolyakGradient(obj,x,~)
            d = obj.designVariable;
            d.update(x);
            obj.cost.computeFunctionAndGradient(d);
        end

        function computeNesterovGradient(obj,~,y)
            d = obj.designVariable;
            d.update(y);
            obj.cost.computeFunctionAndGradient(d);
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

