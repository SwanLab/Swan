classdef OptimizerNullSpace < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'NullSpace';
    end

    properties (Access = private)
        tau
        lineSearchTrials
        lineSearch
        costOld
        upperBound
        lowerBound
        tol = 1e-5
        nX
        hasConverged
        acceptableStep
        oldDesignVariable
        oldCost
        incrementalScheme
        hasFinished
        mOld
        meritNew
        nConstr
        meritGradient
        aJmax
        aGmax
        aJ
        aG
        eta

        globalCost
        globalConstraint
        globalCostGradient
        globalMerit
        globalLineSearch
        globalDual
        globalDesignVar
    end

    methods (Access = public) 
        
        function obj = OptimizerNullSpace(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.outputFunction.monitoring.create(cParams);
            obj.createPrimalUpdater(cParams);
            obj.createDualUpdater(cParams);
            obj.prepareFirstIter();
            obj.aJmax = obj.nullSpaceParameterEstimation(cParams);
            obj.aGmax = obj.rangeSpaceParameterEstimation(cParams);
        end

        function solveProblem(obj)
            obj.hasConverged = false;
            obj.hasFinished = false;
            obj.printOptimizerVariable();
            while ~obj.hasFinished
                obj.update();
                obj.updateIterInfo();
                obj.updateMonitoring();
                obj.checkConvergence();
                obj.printOptimizerVariable();
                obj.checkParameters();
            end
        end

    end

    methods(Access = private)

        function init(obj,cParams)
            obj.upperBound        = cParams.uncOptimizerSettings.ub;
            obj.lowerBound        = cParams.uncOptimizerSettings.lb;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.incrementalScheme = cParams.incrementalScheme;
            obj.nConstr           = cParams.constraint.nSF;
            obj.nX                = obj.designVariable.fun.nDofs;
            obj.maxIter           = cParams.maxIter;
            obj.hasConverged      = false;
            obj.nIter             = 0;
        end

        function aJmax = nullSpaceParameterEstimation(obj,cParams)
            if isfield (cParams.optimizerNames,'aJmax')
                aJmax = cParams.optimizerNames.aJmax;
            else
                DJ = obj.cost.gradient;
                Dg = obj.constraint.gradient;
                aJmax = -1/((Dg'*Dg)\Dg'*DJ);
            end
        end

        function aGmax = rangeSpaceParameterEstimation(obj,cParams)
            if isfield (cParams.optimizerNames,'aGmax')
                aGmax = cParams.optimizerNames.aGmax;
            else
                Dg = obj.constraint.gradient;
                aGmax = 150/(inv(Dg'*Dg));
            end
        end

        function checkParameters(obj)
            if abs(obj.meritNew - obj.mOld) < 10*obj.tol
                g  = obj.constraint.value;
                DJ = obj.cost.gradient;
                if obj.aG <= 0.5*obj.aGmax
                    exponent = -sign(g)*sign(sum(DJ));
                    obj.aJmax = obj.aJmax*1^exponent;
                else
                    exponent  = 1-sign(obj.checkConstraint());
                    obj.aGmax = obj.aGmax*1^exponent;
                end
            end
        end

        function prepareFirstIter(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.costOld = obj.cost.value;
            obj.designVariable.updateOld();
            obj.dualVariable.value = zeros(obj.nConstr,1);
        end

        function updateNullSpaceCoefficient(obj)
            targetVolume = obj.targetParameters.Vfrac;
            obj.aJ       = obj.aJmax*(1-targetVolume);
        end

        function updateRangeSpaceCoefficient(obj)
            obj.aG = obj.aGmax;
        end

        function updateMaximumVolumeRemoved(obj)
            obj.eta = 0.05;
        end

        function update(obj)
            obj.updateNullSpaceCoefficient();
            obj.updateRangeSpaceCoefficient();
            obj.updateMaximumVolumeRemoved();
            x0 = obj.designVariable.fun.fValues;
            g0 = obj.constraint.value;
            obj.saveOldValues(x0);
            obj.calculateInitialStep();
            obj.acceptableStep      = false;
            obj.lineSearchTrials    = 0;
            d.nullSpaceCoefficient  = obj.aJ;
            d.rangeSpaceCoefficient = obj.aG;
            obj.dualUpdater.update(d);
            obj.mOld = obj.computeMeritFunction();
            obj.computeMeritGradient();

            while ~obj.acceptableStep
                x = obj.updatePrimal();
                obj.designVariable.update(x);
                obj.cost.computeFunctionAndGradient();
                obj.constraint.computeFunctionAndGradient();
                s.x0 = x0;
                s.g0 = g0;
                obj.checkStep(s);
            end
        end

        function calculateInitialStep(obj)
            x  = obj.designVariable.fun.fValues;
            DJ = obj.cost.gradient;
            if obj.nIter == 0
                factor = 2e6;
                obj.primalUpdater.computeFirstStepLength(DJ,x,factor);
            else
                factor = 1.2;
                obj.primalUpdater.increaseStepLength(factor);
            end
        end

        function x = updatePrimal(obj)
            x       = obj.designVariable.fun.fValues;
            g       = obj.meritGradient;
            x       = obj.primalUpdater.update(g,x);
        end

        function computeMeritGradient(obj)
            DJ   = obj.cost.gradient;
            Dg   = obj.constraint.gradient;
            l   = obj.dualVariable.value;
            DmF = DJ+Dg*l;
            obj.meritGradient = DmF;
        end

        function checkStep(obj,s)
            x0   = s.x0;
            g0   = s.g0;
            mNew = obj.computeMeritFunction();
            g    = obj.constraint.value;
            v0   = obj.computeVolume(g0);
            v    = obj.computeVolume(g);
            if mNew < obj.mOld && norm(v-v0) < obj.eta
                obj.acceptableStep = true;
                obj.meritNew = mNew;
                obj.dualUpdater.updateOld();
            elseif obj.primalUpdater.isTooSmall()
                warning('Convergence could not be achieved (step length too small)')
                obj.acceptableStep = true;
                obj.meritNew = obj.mOld;
                obj.designVariable.update(x0);
                obj.dualUpdater.updateOld();
                switch obj.designVariable.type
                    case 'Density'
                        obj.primalUpdater.tau = 6000;
                    case 'LevelSet'
                        obj.primalUpdater.tau = 0.01;
                end
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function v = computeVolume(obj,g)
            targetVolume = obj.targetParameters.Vfrac;
            v            = targetVolume*(1+g);
        end

        function mF = computeMeritFunction(obj)
            l  = obj.dualVariable.value;
            J  = obj.cost.value;
            h  = obj.constraint.value;
            mF = J+l'*h;
        end

        function obj = saveOldValues(obj,x)
            obj.oldCost            = obj.cost.value;
            obj.oldDesignVariable  = x;
        end

        function obj = updateOldValues(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
        end

        function obj = checkConvergence(obj)
           if abs(obj.meritNew - obj.mOld) < obj.tol && obj.checkConstraint()
               obj.hasConverged = true;
               if obj.primalUpdater.isTooSmall()
                   obj.primalUpdater.tau = 1;
               end
           end
        end

        function obj = updateMonitoring(obj)
            s.nIter            = obj.nIter;
            s.tau              = obj.primalUpdater.tau;
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
            if iStep < nStep
                itHas = obj.nIter >= iStep;
            else
                itHas = obj.nIter >= obj.maxIter;
            end
        end

    end

end