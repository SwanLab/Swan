classdef OptimizerBisection < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'Bisection';
    end

    properties (Access = private)
        tau
        lineSearchTrials
        lineSearch
        maxLineSearchTrials = 100
        costOld
        tol = 1e-5
        nX
        hasConverged
        acceptableStep
        oldDesignVariable
        oldCost
        problem
        options
        hasFinished
        mOld
        meritNew
        constrProjector
        isInitialStep
        Vtar

        globalCost
        globalConstraint
        globalCostGradient
        globalMerit
        globalLineSearch
        globalDual
    end

    methods (Access = public) 
        
        function obj = OptimizerBisection(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.createPrimalUpdater(cParams);
            s.primalUpdater = obj.primalUpdater;
            obj.constrProjector = ConstraintProjector(cParams,s);
            obj.prepareFirstIter();
        end

        function obj = solveProblem(obj)
            obj.hasConverged  = false;
            obj.isInitialStep = true;
            obj.hasFinished = false;
            obj.printOptimizerVariable();
            d = obj.designVariable;
            obj.constraint.computeFunctionAndGradient(d);
            obj.updateMonitoring();
            while ~obj.hasFinished
                obj.update();
                obj.updateIterInfo();
                obj.printOptimizerVariable();
                obj.updateMonitoring();
                obj.checkConvergence();
            end
        end

    end

    methods(Access = private)

        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable   = cParams.dualVariable;
            obj.nX             = obj.designVariable.fun.nDofs;
            obj.maxIter        = cParams.maxIter;
            obj.nIter          = 0;
            obj.Vtar           = cParams.volumeTarget;
            obj.createMonitoring(cParams);
        end

        function createMonitoring(obj,cParams)
            titlesF       = obj.cost.getTitleFields();
            titlesConst   = obj.constraint.getTitleFields();
            nSFCost       = length(titlesF);
            titles        = [{'Cost'};titlesF;titlesConst;{'Norm L2 x';'\lambda'}];
            titles  = [titles;{'Volume';'Line Search';'Line Search trials'}];
            chCost = cell(1,nSFCost);
            for i = 1:nSFCost
                chCost{i} = 'plot';
            end
            chartTypes = [{'plot'},chCost,{'plot','log','plot','plot','bar','bar'}];
            s.shallDisplay = cParams.monitoring;
            s.maxNColumns  = 5;
            s.titles       = titles;
            s.chartTypes   = chartTypes;
            obj.monitoring = Monitoring(s);
        end

        function updateMonitoring(obj)
            data = obj.cost.value;
            data = [data;obj.cost.getFields(':')];
            data = [data;obj.constraint.value];
            data = [data;obj.designVariable.computeL2normIncrement()];
            data = [data;obj.dualVariable.fun.fValues];
            data = [data;obj.computeVolume(obj.constraint.value)]; % millorar
            if obj.nIter == 0
                data = [data;0;0];
            else
                data = [data;obj.primalUpdater.tau;obj.lineSearchTrials];
            end
            obj.monitoring.update(obj.nIter,data);
        end

        function prepareFirstIter(obj)
            x = obj.designVariable;
            obj.cost.computeFunctionAndGradient(x);
            obj.costOld = obj.cost.value;
            obj.designVariable.updateOld();
            obj.dualVariable.fun.fValues = 0;
        end

        function obj = update(obj)
            x0 = obj.designVariable.fun.fValues;
            obj.saveOldValues(x0);
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            while ~obj.acceptableStep
                obj.constrProjector.project();
                obj.checkStep();
            end
        end

        function obj = calculateInitialStep(obj)
            x = obj.designVariable;
            obj.cost.computeFunctionAndGradient(x);
            obj.constraint.computeFunctionAndGradient(x);
            l       = obj.dualVariable.fun.fValues;
            DJ      = obj.cost.gradient;
            Dg      = obj.constraint.gradient;
            DmF     = DJ + l*Dg;
            if obj.nIter == 0
                factor = 10;
                obj.primalUpdater.computeFirstStepLength(DmF,x,factor);
            else
                obj.tau = 1.5*obj.tau;
            end
        end

        function checkStep(obj)
            x = obj.designVariable;
            obj.cost.computeFunctionAndGradient(x);
            J = obj.cost.value;
            if obj.isInitialStep
                obj.acceptableStep = true;
                obj.isInitialStep  = false;
            else
                if J < obj.oldCost
                    obj.acceptableStep = true;
                    factor = 1.5;
                    obj.primalUpdater.increaseStepLength(factor);
                elseif obj.tau < 1e-10
                    error('Convergence could not be achieved (step length too small)')
                else
                    obj.primalUpdater.decreaseStepLength();
                    obj.lineSearchTrials = obj.lineSearchTrials + 1;
                end
            end
        end

        function v = computeVolume(obj,g)
            targetVolume = obj.Vtar;
            v            = targetVolume*(1+g);
        end

        function mF = computeMeritFunction(obj,x)
            obj.designVariable.update(x)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            J  = obj.cost.value;
            g  = obj.constraint.value;
            l  = obj.lambda;
            mF = J + l*g;
        end

        function obj = saveOldValues(obj,x)
            obj.designVariable.update(x);
            d = obj.designVariable;
            obj.cost.computeFunctionAndGradient(d);
            obj.constraint.computeFunctionAndGradient(d);
            obj.oldCost            = obj.cost.value;
            obj.oldDesignVariable  = x;
        end

        function obj = updateOldValues(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
        end

        function obj = checkConvergence(obj)
           if abs(obj.oldCost - obj.cost.value) < obj.tol
               obj.hasConverged = true;
           else
               
           end

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

        function saveVariablesForAnalysis(obj)
            i                         = obj.nIter + 1;
            obj.globalCost(i)         = obj.cost.value;
            obj.globalConstraint(i)   = obj.constraint.value;
            obj.globalCostGradient(i) = norm(obj.cost.gradient);
            obj.globalMerit(i)        = obj.cost.value;
            obj.globalDual(i)         = obj.dualVariable.fun.fValues;
        end

    end

end