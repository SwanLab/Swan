classdef AcceleratedGradientDescent < handle

    properties (Access = public)
        xFinal
        incJ
        J
        costFields
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
        solverTol
        monitoring
        eta
        matrixFree
    end

    methods (Access = public)

        function obj = AcceleratedGradientDescent(cParams)
            obj.init(cParams);
            obj.createMonitoring();
        end

        function solve(obj)
            if obj.matrixFree
                obj.solveWithMatrixFree();
            else
                obj.solveStandard();
            end
            % obj.solveStandard();
        end

        function solveStandard(obj)
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
                d.updateOld();
                x    = xNew;
                d.update(x);
                iNorm = obj.eta*d.computeNonScaledL2normIncrement();
                obj.solverTol.compute(iNorm);
                obj.cost.computeFunctionAndGradient(d);
                obj.storeIterationData();
                obj.oldCost = obj.cost.value;
                obj.nIter   = obj.nIter + 1;
                obj.updateMonitoring();
                d.plot();
            end
            obj.xFinal = x;
        end

        function solveWithMatrixFree(obj)
            d = obj.designVariable;
            obj.restartDesignVariable();
            x    = obj.designVariable.fun.fValues;
            xOld = x;
            obj.cost.computeFunctionAndGradient(d);
            obj.oldCost = 0;
            obj.nIter   = 1;
            obj.storeIterationData();
            t  = obj.tau;
            t0 = obj.tau;
            while ~obj.hasConverged()
                beta = obj.computeBeta();
                y    = obj.addMomentumTerm(xOld,x,beta);
                obj.gradientComputer(x,y);
                DJ   = obj.cost.gradient;
                xNew = obj.computeGradientStep(y,DJ,t);
                xNew = obj.computeProjection(xNew);
                xOld = x;
                d.updateOld();
                x    = xNew;
                d.update(x);
                % while d.computeNonScaledL2normIncrement() > 1e-4
                %     % d.update(xOld);
                %     t = t/1.5;
                %     xNew = obj.computeGradientStep(y,DJ,t);
                %     xNew = obj.computeProjection(xNew);
                %     d.update(xNew);
                % end
                % t = t0;
                iNorm = obj.eta*d.computeNonScaledL2normIncrement();
                obj.solverTol.compute(iNorm);
                obj.cost.computeFunctionAndGradient(d);
                obj.storeIterationData();
                obj.oldCost = obj.cost.value;
                obj.nIter   = obj.nIter + 1;
                obj.updateMonitoring();
                d.plot();
            end
            obj.xFinal = x;
        end

    end

    methods (Access = private)
        function createMonitoring(obj)
            titles = {'Cost';'Compliance';'Volume';'Norm L2 x'};
            chartTypes = {'plot','plot','plot','log'};
            s.shallDisplay = true;
            s.maxNColumns  = 2;
            s.titles       = titles;
            s.chartTypes   = chartTypes;
            obj.monitoring = Monitoring(s);
        end

        function updateMonitoring(obj)
            data = [obj.cost.value;obj.cost.getFields(1);obj.cost.getFields(2);obj.designVariable.computeNonScaledL2normIncrement()];
            obj.monitoring.update(obj.nIter,data);
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
            obj.solverTol         = cParams.solverTol;
            obj.eta               = cParams.eta;
            obj.momentumParameter = MomentumParameter.create(cParams.momentumParameter);
            obj.gradientComputer  = obj.createGradient(cParams.gDescentType);
            obj.matrixFree        = cParams.matrixFree;
        end

        function storeIterationData(obj)
            obj.incJ(obj.nIter) = abs(obj.cost.value - obj.oldCost);
            obj.J(obj.nIter)    = obj.cost.value;
            obj.costFields(:,obj.nIter) = [obj.cost.getFields(1);obj.cost.getFields(2)];
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

