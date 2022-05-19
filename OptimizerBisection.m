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
        upperBound
        lowerBound
        tol = 1e-3
        nX
        hasConverged
        acceptableStep
        oldDesignVariable
        oldCost
        problem
        options
        incrementalScheme
        hasFinished
        mOld
        meritNew
        constrProjector
        isInitialStep

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
            obj.constrProjector = ConstraintProjector(cParams);
            obj.outputFunction.monitoring.create(cParams);
            obj.prepareFirstIter();
        end

        function obj = solveProblem(obj)
            obj.hasConverged  = false;
            obj.isInitialStep = true;
            while ~obj.hasConverged
                obj.update();
                obj.updateIterInfo();
                obj.updateMonitoring();
                obj.checkConvergence();
                obj.saveVariablesForAnalysis();
            end
        end

    end

    methods(Access = private)

        function init(obj,cParams)
            obj.upperBound             = cParams.uncOptimizerSettings.ub;
            obj.lowerBound             = cParams.uncOptimizerSettings.lb;
            obj.cost                   = cParams.cost;
            obj.constraint             = cParams.constraint;
            obj.designVariable         = cParams.designVar;
            obj.dualVariable           = cParams.dualVariable;
            obj.incrementalScheme      = cParams.incrementalScheme;
            obj.nX                     = length(obj.designVariable.value);
            obj.maxIter                = cParams.maxIter;
            obj.nIter                  = 0;
        end

        function prepareFirstIter(obj)
            obj.cost.computeFunctionAndGradient();
            obj.costOld = obj.cost.value;
            obj.designVariable.updateOld();
            obj.dualVariable.value = 0;
        end

        function obj = update(obj)
            x0 = obj.designVariable.value;
            obj.saveOldValues(x0);
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            while ~obj.acceptableStep
                t = obj.tau;
                obj.constrProjector.project(t);
                obj.checkStep();
            end
        end

        function obj = calculateInitialStep(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            x       = obj.designVariable.value;
            l       = obj.dualVariable.value;
            DJ      = obj.cost.gradient;
            Dg      = obj.constraint.gradient;
            DmF     = DJ + l*Dg;
            if obj.nIter == 0
                obj.tau = 10*sqrt(norm(DmF)/norm(x));
            else
                obj.tau = 1.5*obj.tau;
            end
        end

        function checkStep(obj)
            obj.cost.computeFunctionAndGradient();
            J = obj.cost.value;
            if obj.isInitialStep
                obj.acceptableStep = true;
                obj.isInitialStep  = false;
            else
                if J < obj.oldCost
                    obj.acceptableStep = true;
                elseif obj.tau < 1e-10
                    error('Convergence could not be achieved (step length too small)')
                else
                    obj.tau = obj.tau/2;
                    obj.lineSearchTrials = obj.lineSearchTrials + 1;
                end
            end
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
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.oldCost            = obj.cost.value;
            obj.oldDesignVariable  = x;
        end

        function obj = updateOldValues(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
        end

        function obj = checkConvergence(obj)
           if abs(obj.oldCost - obj.cost.value) < obj.tol && abs(obj.constraint.value) <= 1e-4
               obj.hasConverged = true;
           else
               
           end

        end

        function obj = updateMonitoring(obj)
            s.nIter            = obj.nIter;
            s.tau              = obj.tau;
            s.lineSearch       = obj.lineSearch;
            s.lineSearchTrials = obj.lineSearchTrials;
            s.oldCost          = obj.oldCost;
            s.hasFinished      = obj.hasFinished;
            obj.outputFunction.monitoring.compute(s);
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
            iStep = obj.incrementalScheme.iStep;
            nStep = obj.incrementalScheme.nSteps;
            itHas = obj.nIter >= obj.maxIter*(iStep/nStep);
        end

        function saveVariablesForAnalysis(obj)
            i                         = obj.nIter + 1;
            obj.globalCost(i)         = obj.cost.value;
            obj.globalConstraint(i)   = obj.constraint.value;
            obj.globalCostGradient(i) = norm(obj.cost.gradient);
            obj.globalMerit(i)        = obj.cost.value;
            obj.globalLineSearch(i)   = obj.tau;
            obj.globalDual(i)         = obj.dualVariable.value;
            iStep                     = obj.incrementalScheme.iStep;
            nStep                     = obj.incrementalScheme.nSteps;
            if obj.hasConverged && iStep == nStep
                c = obj.globalCost;
                h = obj.globalConstraint;
                g = obj.globalCostGradient;
                m = obj.globalMerit;
                t = obj.globalLineSearch;
                d = obj.globalDual;
                save('BisectionVariables.mat',"t","m","c","g","h","d");
            end
        end

    end

end