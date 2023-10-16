classdef OptimizerInteriorPoint < Optimizer

    properties (Access = private)
        slack
        bounds
        pusher
        tol = 1e-5
        barrierTau
        lineSearch
        lineSearchTrials
        error
        hessian
        hessianUpdated
        dx
        ds 
        dlam
        xNew
        sNew
        hasConverged
        hasFinished
        incrementalScheme
        nSlack
        baseVariables
        oldDesignVariable
        oldCost
        oldCostGradient
        mOld
        meritNew
        acceptableStep
    end

    methods (Access = public)

        function obj = OptimizerInteriorPoint(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
            obj.outputFunction.monitoring.create(cParams);
            obj.createPrimalUpdater(cParams);
            obj.createDualUpdater(cParams);
        end

        function solveProblem(obj)
            obj.hasConverged = false;
            obj.hasFinished = false;
            obj.previousComputations();
            while ~obj.hasFinished
                obj.update();
                obj.updateIterInfo();
                obj.updateMonitoring();
                obj.checkConvergence();
                obj.checkNewBarrierProblem();
            end
        end
    end

    methods(Access = private)
        function createVariablesPusher(obj)
            s.tau            = obj.barrierTau;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.slack          = obj.slack;
            s.bounds         = obj.bounds;
            obj.pusher       = IPMVariablesPusher(s);
        end

        function previousComputations(obj)
            obj.computeInitialBarrierParameter();
            obj.createVariablesPusher();
            obj.moveVariablesToFeasible();
            obj.computeInitialDualVariableBounds();
            obj.updateCostAndConstraintWithPenalty();
            obj.dualUpdater.compute(obj.bounds.zLB,obj.bounds.zUB);
        end

        function moveVariablesToFeasible(obj)
            x         = obj.pusher.replaceOutOfDesignVarBounds();
            obj.designVariable.update(x);
            obj.slack = obj.pusher.replaceOutOfSlackBounds();
        end

        function computeInitialDualVariableBounds(obj)
            s.mu             = obj.baseVariables.mu;
            s.designVariable = obj.designVariable;
            s.slack          = obj.slack;
            s.bounds         = obj.bounds;
            obj.bounds.zLB   = obj.dualUpdater.computeLowerBound(s);
            obj.bounds.zUB   = obj.dualUpdater.computeUpperBound(s);
        end

        function computeInitialBarrierParameter(obj)
            mu             = obj.baseVariables.mu;
            tauMax         = obj.baseVariables.tau_max;
            obj.barrierTau = min(tauMax,100*mu);
        end

        function updateCost(obj)
            obj.penalizeDueToDesignVariable();
            obj.penalizeDueToSlackVariable();
        end

        function penalizeDueToDesignVariable(obj)
            c              = obj.cost.value;
            mu             = obj.baseVariables.mu;
            s.field        = obj.designVariable.value;
            s.lowerBound   = obj.bounds.xLB';
            s.upperBound   = obj.bounds.xUB';
            obj.cost.value = c - mu*obj.computeLogPenaltyTerm(s);
        end

        function penalizeDueToSlackVariable(obj)
            c              = obj.cost.value;
            mu             = obj.baseVariables.mu;
            s.field        = obj.slack;
            s.lowerBound   = obj.bounds.constraintLB(obj.isInequality());
            s.upperBound   = obj.bounds.constraintUB(obj.isInequality());
            obj.cost.value = c - mu*obj.computeLogPenaltyTerm(s);
        end

        function itIs = isInequality(obj)
            itIs = obj.bounds.constraintUB > obj.bounds.constraintLB;
        end

        function updateCostGradient(obj)
            nS                = size(obj.slack,2);
            obj.cost.gradient = [obj.cost.gradient' zeros(1,nS)];
        end

        function updateConstraint(obj)
            g                      = obj.constraint.value;
            gLB                    = obj.bounds.constraintLB(:,1);
            s                      = obj.slack;
            g(~obj.isInequality()) = g(~obj.isInequality()) - gLB(~obj.isInequality());
            g(obj.isInequality())  = g(obj.isInequality()) - s;
            obj.constraint.value   = g;
        end

        function updateConstraintGradient(obj)
            ndof    = size(obj.constraint.gradient,1);
            iSl     = ndof+1:ndof+obj.nSlack;
            iConstr = obj.isInequality()==true;
            obj.constraint.gradient(iSl,iConstr) = -1;
        end

        function computeHessian(obj)
            s.initHessian             = obj.hessian;
            s.designVariable.value    = obj.designVariable.value;
            s.designVariable.valueOld = obj.oldDesignVariable;
            s.cost.gradient           = obj.cost.gradient;
            s.cost.gradientOld        = obj.oldCostGradient;
            h                         = HessianComputer(s);
            obj.hessian               = h.compute();
        end

        function meritF = computeMeritFunction(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            J      = obj.cost.value;
            g      = obj.constraint.value;
            nu     = obj.baseVariables.nu;
            meritF = J + nu*sum(abs(g));
            obj.updateCostAndConstraintWithPenalty();
        end

        function updateCostAndConstraintWithPenalty(obj)
            obj.updateCost();
            obj.updateCostGradient();
            obj.updateConstraint();
            obj.updateConstraintGradient();
        end

        function update(obj)
            x0 = obj.designVariable.value;
            obj.saveOldValues(x0);
            obj.mOld = obj.computeMeritFunction(x0);
            obj.computeOptimizerDirections();
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;            
            while ~obj.acceptableStep
                obj.xNew = obj.updatePrimal();
                obj.updatePrimalVariables();
                obj.updateDual();
                obj.pushVariables();
                obj.checkStep(obj.xNew,x0);
            end
        end

        function obj = calculateInitialStep(obj)
            x  = obj.designVariable.value;
            DJ = obj.cost.gradient;
            if obj.nIter == 0
                factor = 1;
                obj.primalUpdater.computeFirstStepLength(DJ,x,factor);
            else
                factor = 1.05;
                obj.primalUpdater.increaseStepLength(factor);
            end
        end

        function x = updatePrimal(obj)
            x = obj.designVariable.value;
            g = -obj.dx;
            x = obj.primalUpdater.update(g,x);
        end

        function updatePrimalVariables(obj)
            alphaPrimal  = obj.primalUpdater.tau;
            if obj.nSlack >= 1
                obj.sNew = obj.slack + alphaPrimal * obj.ds';
            else
                obj.sNew = [];
            end
        end

        function updateDual(obj)
            obj.dualUpdater.updateAlpha(obj.primalUpdater.tau);
            obj.dualUpdater.update(obj.dlam);
            obj.bounds.zNewLB = obj.dualUpdater.updateLowerBound(obj.bounds);
            obj.bounds.zNewUB = obj.dualUpdater.updateLowerBound(obj.bounds);  
            obj.pusher.updateBounds(obj.bounds);
        end

        function pushVariables(obj)
            obj.xNew   = obj.pusher.pushDesignVariable(obj.xNew);
            obj.sNew   = obj.pusher.pushSlack();
            obj.bounds = obj.pusher.pushDualVariableBounds();
        end

        function updateWithAcceptance(obj)
            obj.designVariable.update(obj.xNew);
            obj.slack      = obj.sNew;
            obj.bounds.zLB = obj.bounds.zNewLB;
            obj.bounds.zUB = obj.bounds.zNewUB;
        end

        function obj = saveOldValues(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.oldCost             = obj.cost.value;
            obj.oldCostGradient     = obj.cost.gradient;
            obj.designVariable.updateOld;
            obj.oldDesignVariable   = obj.designVariable.value;
        end

        function checkStep(obj,x,x0)
            mNew                      = obj.computeMeritFunction(x);
            obj.updateNuValue();
            if mNew < obj.mOld
                obj.acceptableStep    = true;
                obj.updateWithAcceptance();
                obj.meritNew          = mNew;
            elseif obj.primalUpdater.isTooSmall()
                warning('Convergence could not be achieved (step length too small)')
                obj.acceptableStep    = true;
                obj.primalUpdater.tau = 1;
                obj.meritNew          = mNew;
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials  = obj.lineSearchTrials + 1;
            end
        end

        function updateNuValue(obj)
            DJ                     = obj.cost.gradient;
            Dx                     = [obj.dx;obj.ds];
            hess                   = obj.hessianUpdated;
            costDecreasePrediction = DJ*Dx + max(0,0.5*Dx'*hess*Dx);
            theta                  = sum(abs(obj.constraint.value));
            r                      = 0.1;
            nuUpdated              = costDecreasePrediction/((1-r)*theta);
            obj.baseVariables.nu   = max(1,min(1000,nuUpdated));
        end

        function computeOptimizerDirections(obj)
            s.cost             = obj.cost;
            s.constraint       = obj.constraint;
            s.designVariable   = obj.designVariable;
            s.slack            = obj.slack;
            s.dualVariable     = obj.dualVariable;
            s.baseVariables    = obj.baseVariables;
            s.hessian          = obj.hessian;
            s.nSlack           = obj.nSlack;
            s.bounds           = obj.bounds;
            dirs               = IPMDirectionComputer(s);
            dirs.compute();
            obj.dx             = dirs.gradients.dx;
            obj.ds             = dirs.gradients.ds;
            obj.dlam           = dirs.gradients.dlam;
            obj.bounds.dzLB    = dirs.gradients.dzL;
            obj.bounds.dzUB    = dirs.gradients.dzU;
            obj.hessianUpdated = dirs.updatedHessian;
        end
        
        function init(obj,cParams)
            obj.bounds.xUB        = cParams.uncOptimizerSettings.ub;
            obj.bounds.xLB        = cParams.uncOptimizerSettings.lb;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.incrementalScheme = cParams.incrementalScheme;
            obj.maxIter           = cParams.maxIter;
            obj.nIter             = 0;
            obj.constraintCase    = cParams.constraintCase;
            obj.cost.computeFunctionAndGradient();
            obj.hessian = eye(obj.designVariable.mesh.nnodes);
            obj.constraint.computeFunctionAndGradient();
            obj.loadIPMVariables();
        end

        function loadIPMVariables(obj)
            obj.baseVariables.nu         = 1e6;
            obj.baseVariables.mu         = 1e-9;
            obj.baseVariables.slack_init = true;
            obj.baseVariables.tau_max    = 0.1;
            obj.baseVariables.k_mu       = 0.2;
            obj.computeLinearBounds();
            obj.computeSlackVariables();
            obj.correctUpperAndLowerBounds();
        end

        function computeLinearBounds(obj)
            for i = 1:length(obj.constraintCase)
                if strcmp(obj.constraintCase{i},'INEQUALITY')
                    obj.bounds.constraintLB(i,:) = -inf;
                    obj.bounds.constraintUB(i,:) = 0;
                else
                    obj.bounds.constraintLB(i,:) = 0;
                    obj.bounds.constraintUB(i,:) = 0;
                end
            end
        end

        function correctUpperAndLowerBounds(obj)
            nnode          = obj.designVariable.mesh.nnodes;
            obj.bounds.xLB = (obj.bounds.xLB-1e-6)*ones(1,nnode);
            obj.bounds.xUB = (obj.bounds.xUB+1e-6)*ones(1,nnode);
        end

        function computeSlackVariables(obj)
            obj.bounds.sLB = obj.bounds.constraintLB(obj.isInequality());
            obj.bounds.sUB = obj.bounds.constraintUB(obj.isInequality());
            obj.slack      = zeros(size(obj.bounds.sLB));
            obj.nSlack     = length(obj.slack);
            obj.updateConstraint();
            if (obj.baseVariables.slack_init)
                obj.slack = obj.constraint.value(obj.isInequality());
            else
                obj.slack(:) = 0.01;
            end
        end

        function checkConvergence(obj)
            obj.cost.computeFunctionAndGradient();
            obj.computeHessian();
            obj.constraint.computeFunctionAndGradient();
            obj.updateCostAndConstraintWithPenalty();
            obj.computeError();
            if obj.error <= obj.tol
                obj.hasConverged = true;
            end
        end

        function computeError(obj)
            s.sMax           = 100;
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.dualVariable   = obj.dualVariable;
            s.slack          = obj.slack;
            s.bounds         = obj.bounds;
            eps              = IPMErrorComputer(s);
            eps.compute();
            obj.error        = eps.error;
        end

        function checkNewBarrierProblem(obj)
            k_mu = obj.baseVariables.k_mu;
            if (obj.error < k_mu * obj.baseVariables.mu)
                obj.baseVariables.mu = k_mu*obj.baseVariables.mu;
                obj.barrierTau       = min(obj.baseVariables.tau_max,100*obj.baseVariables.mu);
            end
        end

        function obj = updateMonitoring(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            s.nIter            = obj.nIter;
            s.tau              = obj.primalUpdater.tau;
            s.lineSearch       = obj.lineSearch;
            s.lineSearchTrials = obj.lineSearchTrials;
            s.oldCost          = obj.oldCost;
            s.hasFinished      = obj.hasFinished;
            s.meritNew         = obj.meritNew;
            s.lambda           = obj.dualVariable.value;
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

    methods (Access = private,Static)

        function penalty = computeLogPenaltyTerm(s)
            x       = s.field;
            lb      = s.lowerBound;
            ub      = s.upperBound;
            penalty = sum(log(x-lb) + log(ub-x));
        end

    end

end