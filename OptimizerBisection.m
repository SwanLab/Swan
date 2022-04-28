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
            while ~obj.hasConverged
                obj.update();
                obj.updateIterInfo();
                obj.updateMonitoring();
                obj.checkConvergence();
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
            obj.hasConverged           = false;
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
            obj.mOld = obj.computeMeritFunction(x0);
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            while ~obj.acceptableStep
                x = obj.updatePrimal();
                x = obj.constrProjector.project(x);
                obj.checkStep(x,x0);
            end
            obj.updateOldValues(x);
        end

        function obj = calculateInitialStep(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            x       = obj.designVariable.value;
            l       = obj.lambda;
            DJ      = obj.cost.gradient;
            Dg      = obj.constraint.gradient;
            DmF     = DJ + l*Dg;
            if obj.nIter == 0
                obj.tau = 1*sqrt(norm(DmF)/norm(x));
            else
                obj.tau = 1.5*obj.tau;
            end
        end

        function obj = updateDualDirect(obj)
            obj.constraint.computeFunctionAndGradient();
            obj.cost.computeFunctionAndGradient();
            DJ = obj.cost.gradient;
            Dg = obj.constraint.gradient;
            g  = obj.constraint.value;
            S  = (Dg'*Dg)^-1;
            t  = obj.tau;
            aC = 1;
            aJ = 1;
            l  = aC/aJ*S*(g - t*Dg'*DJ);
            obj.lambda = l;
        end

        function x = updatePrimal(obj)
            lb      = obj.lowerBound;
            ub      = obj.upperBound;
            t       = obj.tau;
            Dg      = obj.constraint.gradient;
            DJ      = obj.cost.gradient;
            l       = obj.lambda;
            x       = obj.designVariable.value;
            dx      = -t*(DJ + l*Dg);
            xN      = x + dx;
            x       = min(ub,max(xN,lb));
        end

        function checkStep(obj,x,x0)
            mNew = obj.computeMeritFunction(x);

            if obj.nIter == 0 && mNew == obj.mOld
                obj.acceptableStep = true;
                obj.updateDualDirect();
                obj.meritNew = mNew;
            end
            if mNew < obj.mOld
                obj.acceptableStep = true;
                obj.updateDualDirect();
                obj.meritNew = mNew;
            elseif obj.tau < 1e-10
                error('Convergence could not be achieved (step length too small)')
            else
                obj.tau = obj.tau/2;
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
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
           if abs(obj.oldCost - obj.cost.value) < obj.tol && max(obj.constraint.value) <= 0
               obj.hasConverged = true;
           else
               
           end

        end

        function obj = updateMonitoring(obj)
            obj.updateIterInfo();
            s.nIter            = obj.nIter;
            s.tau              = obj.tau;
            s.lineSearch       = obj.lineSearch;
            s.lineSearchTrials = obj.lineSearchTrials;
            s.oldCost          = obj.oldCost;
            s.hasFinished      = obj.hasFinished;
            s.meritNew         = obj.meritNew;
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

    end

end