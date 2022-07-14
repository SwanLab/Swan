classdef OptimizerAugmentedLagrangian < Optimizer

    properties (GetAccess = public, SetAccess = protected)
        type = 'Augmented Lagrangian';
    end

    properties (Access = private)
        tau
        lineSearchTrials
        lineSearch
        costOld
        upperBound
        lowerBound
        tol = 1e-8
        nX
        nConstr
        hasConverged
        acceptableStep
        oldDesignVariable
        oldCost
        incrementalScheme
        hasFinished
        mOld
        meritNew
        penalty
        meritGradient

        globalCost
        globalConstraint
        globalCostGradient
        globalMerit
        globalLineSearch
        globalDual
        globalDesignVar
    end

    methods (Access = public) 
        
        function obj = OptimizerAugmentedLagrangian(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.outputFunction.monitoring.create(cParams);
            obj.createPrimalUpdater(cParams);
            obj.createDualUpdater(cParams);
            obj.prepareFirstIter();
        end

        function obj = solveProblem(obj)
            obj.hasConverged = false;
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
%             obj.saveVariablesForAnalysis();
            while ~obj.hasConverged
                obj.update();
                obj.updateIterInfo();
                obj.updateMonitoring();
                obj.checkConvergence();
%                 obj.saveVariablesForAnalysis();
            end
        end

    end

    methods(Access = private)

        function init(obj,cParams)
            obj.upperBound             = cParams.uncOptimizerSettings.ub;
            obj.lowerBound             = cParams.uncOptimizerSettings.lb;
            obj.cost                   = cParams.cost;
            obj.constraint             = cParams.constraint;
            obj.nConstr                = cParams.constraint.nSF;
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
            obj.dualVariable.value = zeros(obj.nConstr,1);
            obj.penalty            = 10;
        end

        function obj = update(obj)
            x0 = obj.designVariable.value;
            obj.designVariable.update(x0);
            obj.saveOldValues(x0);
            obj.mOld = obj.computeMeritFunction(x0);
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            obj.computeMeritGradient();
            while ~obj.acceptableStep
                x = obj.updatePrimal();
                obj.checkStep(x,x0);
            end
            obj.updateOldValues(x);
        end

        function displayIter(obj,x)
            m = obj.designVariable.mesh.innerMeshOLD;
            bm = m.createBoundaryMesh();
            s.backgroundMesh = m;
            s.boundaryMesh   = bm;
            um = UnfittedMesh(s);
            um.compute(x);
            figure()
            um.plot();
        end

        function obj = calculateInitialStep(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            x       = obj.designVariable.value;
            l       = obj.dualVariable.value;
            DJ      = obj.cost.gradient;
            Dg      = obj.constraint.gradient;
            g       = obj.constraint.value;
            p       = obj.penalty;
            DmF     = DJ + Dg*(l + p*g);
            if obj.nIter == 0
                factor = 1;
                obj.primalUpdater.computeFirstStepLength(DmF,x,factor);
            else
                factor = 1.05;
                obj.primalUpdater.increaseStepLength(factor);
            end
        end

        function x = updatePrimal(obj)
            x   = obj.designVariable.value;
            g   = obj.meritGradient;
            x   = obj.primalUpdater.update(g,x);
        end

        function computeMeritGradient(obj)
            Dh    = obj.constraint.gradient;
            DJ    = obj.cost.gradient;
            l     = obj.dualVariable.value;
            p     = obj.penalty;
            gPlus = obj.defineConstraintValue();
            g     = (DJ + Dh*(l + p*gPlus));
            obj.meritGradient = g;
        end

        function mF = computeMeritFunction(obj,x)
            obj.designVariable.update(x)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            J      = obj.cost.value;
            gPlus  = obj.defineConstraintValue();
            l      = obj.dualVariable.value;
            rho    = obj.penalty;
            mF     = J + l'*gPlus + 0.5*rho*(gPlus'*gPlus);
        end

        function c = defineConstraintValue(obj)
            c   = obj.constraint.value;
            l   = obj.dualVariable.value;
            rho = obj.penalty;
            for i = 1:obj.nConstr
                switch obj.constraintCase{i}
                    case 'EQUALITY'
                        
                    case 'INEQUALITY'
                        c(i) = max(c(i),-l/rho);
                end
            end
        end

        function checkStep(obj,x,x0)
            mNew = obj.computeMeritFunction(x);
            if mNew < obj.mOld
                obj.acceptableStep = true;
                obj.dualUpdater.updatePenalty(obj.penalty);
                obj.dualUpdater.update();
                obj.meritNew = mNew;
            elseif obj.primalUpdater.isTooSmall()
                error('Convergence could not be achieved (step length too small)')
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
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
           if abs(obj.meritNew - obj.mOld) < obj.tol && obj.checkConstraint()
               obj.hasConverged = true;
           else
               
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
            itHas = obj.nIter >= obj.maxIter*(iStep/nStep);
        end

        function saveVariablesForAnalysis(obj)
            i                           = obj.nIter + 1;
            obj.globalCost(i)           = obj.cost.value;
            obj.globalConstraint(:,i)   = obj.constraint.value;
            obj.globalCostGradient(i)   = norm(obj.cost.gradient);
            obj.globalMerit(i)          = obj.meritNew;
            obj.globalLineSearch(i)     = obj.primalUpdater.tau;
            obj.globalDual(:,i)         = obj.dualVariable.value;
            obj.globalDesignVar(:,i)    = obj.designVariable.value;
            if obj.hasConverged
                c = obj.globalCost;
                h = obj.globalConstraint;
                g = obj.globalCostGradient;
                m = obj.globalMerit;
                t = obj.globalLineSearch;
                d = obj.globalDual;
                v = obj.globalDesignVar;
                save('name.mat',"c","g","h","d","v");
            end
        end


    end

end