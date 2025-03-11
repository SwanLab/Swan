classdef OptimizerProjectedGradient < handle

    properties (Access = private)
        tolCost   = 1e-6
    end

    properties (Access = private)
        cost
        designVariable
        primalUpdater
        maxIter
        nIter
        monitoring
        lineSearchTrials
        hasConverged
        acceptableStep
        hasFinished
        meritNew
        meritOld
        meritGradient
        etaNorm
    end

    methods (Access = public) 
        function obj = OptimizerProjectedGradient(cParams)
            obj.init(cParams);
            obj.createMonitoring(cParams);
            obj.prepareFirstIter();
        end

        function solveProblem(obj)
            obj.hasConverged = false;
            obj.hasFinished  = false;
            obj.plotVariable();
            obj.updateMonitoring();
            while ~obj.hasFinished
                obj.update();
                obj.updateIterInfo();
                obj.plotVariable();
                obj.updateMonitoring();
                obj.checkConvergence();
                obj.designVariable.updateOld();
            end
        end
    end

    methods(Access = private)
        function init(obj,cParams)
            obj.cost            = cParams.cost;
            obj.designVariable  = cParams.designVariable;
            obj.maxIter         = cParams.maxIter;
            obj.hasConverged    = false;
            obj.nIter           = 0;
            obj.etaNorm         = cParams.etaNorm;
            obj.createPrimalUpdater(cParams);
        end

        function createPrimalUpdater(obj,cParams)
            s.ub              = cParams.ub;
            s.lb              = cParams.lb;
            s.tauMax          = 100;
            obj.primalUpdater = ProjectedGradient(s);
        end

        function createMonitoring(obj,cParams)
            s.shallDisplay   = cParams.monitoring;
            s.cost           = obj.cost;
            s.designVariable = obj.designVariable;
            s.primalUpdater  = obj.primalUpdater;
            obj.monitoring   = MonitoringProjectedGradient(s);
        end

        function updateMonitoring(obj)
            s.lineSearchTrials = obj.lineSearchTrials;
            s.meritNew         = obj.meritNew;
            obj.monitoring.update(obj.nIter,s);
            obj.monitoring.refresh();
        end

        function plotVariable(obj)
            if ismethod(obj.designVariable,'plot')
                obj.designVariable.plot();
            end
        end

        function prepareFirstIter(obj)
            obj.meritOld = obj.computeMeritFunction();
            obj.designVariable.updateOld();
        end

        function update(obj)
            x0 = obj.designVariable.fun.fValues;
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            obj.computeMeritGradient();
            obj.calculateInitialStep();
            while ~obj.acceptableStep
                obj.updatePrimal();
                obj.checkStep(x0);
            end
        end

        function calculateInitialStep(obj)
            x   = obj.designVariable;
            DmF = obj.meritGradient;
            if obj.nIter == 0
                factor = 50;
                obj.primalUpdater.computeFirstStepLength(DmF,x,factor);
            else
                factor = 1.05;
                obj.primalUpdater.increaseStepLength(factor);
            end
        end

        function updatePrimal(obj)
            x = obj.designVariable;
            g = obj.meritGradient;
            x = obj.primalUpdater.update(g,x);
            obj.designVariable = x;
        end

        function computeMeritGradient(obj)
            DJ = obj.cost.gradient;
            obj.meritGradient = DJ;
        end

        function checkStep(obj,x0)
            mNew = obj.computeMeritFunction();
            x    = obj.designVariable.fun.fValues;
            if mNew <= obj.meritOld+1e-3  &&  norm(x-x0)/(norm(x0)+1) < obj.etaNorm
                obj.acceptableStep = true;
                obj.meritNew       = mNew;
            elseif obj.primalUpdater.isTooSmall()
                warning('Convergence could not be achieved (step length too small)')
                obj.acceptableStep = true;
                obj.meritNew       = obj.mOld;
                obj.designVariable.update(x0);
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function mF = computeMeritFunction(obj)
            x = obj.designVariable;
            obj.cost.computeFunctionAndGradient(x);
            J  = obj.cost.value;
            mF = J;
        end

        function obj = checkConvergence(obj)
            if abs(obj.meritNew - obj.meritOld) < obj.tolCost
                obj.hasConverged = true;
                if obj.primalUpdater.isTooSmall()
                    obj.primalUpdater.tau = 1;
                end
            end
            obj.meritOld = obj.meritNew;
        end

        function updateIterInfo(obj)
            obj.increaseIter();
            obj.updateStatus();
        end

        function increaseIter(obj)
            obj.nIter = obj.nIter + 1;
        end

        function updateStatus(obj)
            obj.hasFinished = obj.hasConverged || obj.hasExceededStepIterations();
        end

        function itHas = hasExceededStepIterations(obj)
            itHas = obj.nIter >= obj.maxIter;
        end
    end
end