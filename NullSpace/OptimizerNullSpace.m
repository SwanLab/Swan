classdef OptimizerNullSpace < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'NullSpace';
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
        nConstr
    end

    methods (Access = public) 
        
        function obj = OptimizerNullSpace(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.outputFunction.monitoring.create(cParams);
            cParams.tau = obj.tau;
            obj.createDualUpdater(cParams);
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
            obj.nConstr                = cParams.nConstr;
            obj.nX                     = length(obj.designVariable.value);
            obj.maxIter                = cParams.maxIter;
            obj.hasConverged           = false;
            obj.nIter                  = 0;
        end

        function prepareFirstIter(obj)
            obj.cost.computeFunctionAndGradient();
            obj.costOld = obj.cost.value;
            obj.designVariable.updateOld();
            obj.dualVariable.value = zeros(obj.nConstr,1);
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
                obj.checkStep(x,x0);
            end
            obj.updateOldValues(x);
        end

        function obj = calculateInitialStep(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            x       = obj.designVariable.value;
            l       = obj.dualVariable.value;
            DJ      = obj.cost.gradient;
            Dg      = obj.constraint.gradient;
            aJ      = 1;
            DmF     = aJ*(DJ + l'*Dg);
            if obj.nIter == 0
                obj.tau = 3*sqrt(norm(DmF)/norm(x));
            else
                obj.tau = 1.05*obj.tau;
            end
            obj.tau = min(obj.tau,1.2);
        end

        function x = updatePrimal(obj)
            lb      = obj.lowerBound;
            ub      = obj.upperBound;
            t       = obj.tau;
            Dg      = obj.constraint.gradient;
            DJ      = obj.cost.gradient;
            l       = obj.dualVariable.value;
            x       = obj.designVariable.value;
            aJ      = 1;
            dAJ     = aJ*(DJ + l'*Dg);
            dx      = -t*(dAJ);
            xN      = x + dx;
            x       = min(ub,max(xN,lb));
        end

        function checkStep(obj,x,x0)
            mNew = obj.computeMeritFunction(x);
            if obj.nIter == 0 && mNew == obj.mOld
                obj.acceptableStep = true;
                obj.dualUpdater.updateTau(obj.tau);
                obj.dualUpdater.update();
                obj.meritNew = mNew;
            end
            if mNew < obj.mOld
                obj.acceptableStep = true;
                obj.dualUpdater.updateTau(obj.tau);
                obj.dualUpdater.update();
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
            J      = obj.cost.value;
            g      = obj.constraint.value;
            l      = obj.dualVariable.value;
            aJ     = 1;
            AJ     = aJ*(J + l'*g);
            mF     = AJ;
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
           if abs(obj.oldCost - obj.cost.value) < obj.tol && obj.checkConstraint()
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