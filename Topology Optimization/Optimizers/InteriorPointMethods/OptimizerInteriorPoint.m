classdef OptimizerInteriorPoint < Optimizer

    properties (Access = private)
        slack
        bounds
        tol = 1e-5
        barrierTau
        lineSearch
        lineSearchTrials
        error
        H
        hessian
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
        function previousComputations(obj)
            obj.moveVariablesToFeasible();
            obj.computeInitialBarrierParameter();
            obj.computeDualVariableBounds();
            obj.updateCostAndConstraintWithPenalty();
            obj.dualUpdater.compute(obj.bounds.zLB,obj.bounds.zUB);
        end

        function moveVariablesToFeasible(obj)
            x         = obj.replaceOutOfDesignVarBounds();
            obj.designVariable.update(x);
            obj.slack = obj.replaceOutOfSlackBounds();
        end

        function computeDualVariableBounds(obj)
            mu = obj.baseVariables.mu;
            x  = obj.designVariable.value';
            lx = obj.bounds.xLB;
            ux = obj.bounds.xUB;
            s  = obj.slack';
            ls = obj.bounds.sLB;
            us = obj.bounds.sUB;
            zL = mu./(x-lx);
            zL = [zL,mu./(s-ls)];
            zU = mu./(ux-x);
            zU = [zU,mu./(us-s)];
            obj.bounds.zLB = zL;
            obj.bounds.zUB = zU;
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
            nS = size(obj.slack,2);
            obj.cost.gradient = [obj.cost.gradient' zeros(1,nS)];
        end

        function updateConstraint(obj)
            g                      = obj.constraint.value;
            gLB                    = obj.bounds.constraintLB(:,1);
            sl                     = obj.slack;
            g(~obj.isInequality()) = g(~obj.isInequality()) - gLB(~obj.isInequality());
            g(obj.isInequality())  = g(obj.isInequality()) - sl;
            obj.constraint.value   = g;
        end

        function updateConstraintGradient(obj)
            mC = size(obj.constraint.gradient,1);
            nC = size(obj.constraint.gradient,2);
            dk = length(find(obj.isInequality()==true));
            i  = mC+1:mC+dk;
            obj.constraint.gradient(i,nC) = -1;
        end

        function computeHessian(obj)
            nnode     = obj.designVariable.mesh.nnodes;
            hess      = obj.hessian; 
            deltaX    = obj.designVariable.value - obj.oldDesignVariable;
            deltaCost = obj.cost.gradient - obj.oldCostGradient;
            if deltaCost == zeros(size(deltaCost))
                obj.hessian = zeros(nnode,nnode);
            else
                costFrac    = (deltaCost*deltaCost')/(deltaCost'*deltaX);
                xFrac       = (hess*deltaX*(hess*deltaX)')/(deltaX'*hess*deltaX);
                obj.hessian = hess + costFrac - xFrac;
            end
        end

        function meritF = computeMeritFunction(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            J      = obj.cost.value;
            g      = obj.constraint.value;
            meritF = J + obj.baseVariables.nu*sum(abs(g));
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
            alphaDual = obj.primalUpdater.tau;
            obj.dualUpdater.update(alphaDual,obj.dlam);
            obj.bounds.zNewLB = obj.bounds.zLB + alphaDual*obj.bounds.dzLB';
            obj.bounds.zNewUB = obj.bounds.zUB + alphaDual*obj.bounds.dzUB';           
        end

        function pushVariables(obj)
            s.field                    = obj.xNew;
            s.lowerBound               = obj.bounds.xLB;
            s.upperBound               = obj.bounds.xUB;
            s.tau                      = obj.barrierTau;
            obj.xNew                   = obj.pushVarsFunction(s);
            s.field                    = obj.slack;
            s.lowerBound               = obj.bounds.sLB;
            s.upperBound               = obj.bounds.sUB;
            obj.sNew                   = obj.pushVarsFunction(s);
            isZLNeg                    = obj.bounds.zNewLB<0;
            isZUNeg                    = obj.bounds.zNewUB<0;
            tau                        = obj.barrierTau;
            obj.bounds.zNewLB(isZLNeg) = tau*(obj.bounds.zLB(isZLNeg));
            obj.bounds.zNewUB(isZUNeg) = tau*(obj.bounds.zUB(isZUNeg));
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
            gradDesVar             = [obj.dx;obj.ds];
            hess                   = obj.H;
            costPrediction(1)      = DJ*gradDesVar;
            costPrediction(2)      = max(0,0.5*gradDesVar'*hess*gradDesVar);
            costDecreasePrediction = sum(costPrediction);
            theta                  = sum(abs(obj.constraint.value));
            rho                    = 0.1;
            nuUpdated              = costDecreasePrediction/((1-rho)*theta);
            obj.baseVariables.nu   = max(1,min(1000,nuUpdated));
        end

        function computeOptimizerDirections(obj)
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.slack          = obj.slack;
            s.dualVariable   = obj.dualVariable;
            s.baseVariables  = obj.baseVariables;
            s.hessian        = obj.hessian;
            s.nSlack         = obj.nSlack;
            s.bounds         = obj.bounds;
            dirs             = IPMDirectionComputer(s);
            dirs.compute();
            obj.dx           = dirs.gradients.dx;
            obj.ds           = dirs.gradients.ds;
            obj.dlam         = dirs.gradients.dlam;
            obj.bounds.dzLB  = dirs.gradients.dzL;
            obj.bounds.dzUB  = dirs.gradients.dzU;
            obj.H            = dirs.updatedHessian;
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

        function x = replaceOutOfDesignVarBounds(obj)
            s.x   = obj.designVariable.value;
            s.xUB = obj.bounds.xUB;
            s.xLB = obj.bounds.xLB;
            x     = obj.computeOutOfBounds(s);
        end

        function x = replaceOutOfSlackBounds(obj)
            s.x   = obj.slack;
            s.xUB = obj.bounds.sUB;
            s.xLB = obj.bounds.sLB;
            x     = obj.computeOutOfBounds(s);
        end

    end

    methods (Access = private,Static)

        function x = computeOutOfBounds(s)
            x          = s.x;
            isLower    = x<=s.xLB;
            isUpper    = x>=s.xUB;
            x(isLower) = min(s.xUB(isLower),s.xLB(isLower)+1e-2);
            x(isUpper) = max(s.xLB(isUpper),s.xUB(isUpper)-1e-2);
        end

        function penalty = computeLogPenaltyTerm(s)
            x       = s.field;
            lb      = s.lowerBound;
            ub      = s.upperBound;
            penalty = sum(log(x-lb) + log(ub-x));
        end

        function x = pushVarsFunction(s)
            tau        = s.tau;
            x          = s.field;
            lb         = s.lowerBound;
            ub         = s.upperBound;
            isLower    = x<lb;
            isUpper    = x>ub;
            dxl        = tau*(x(isLower) - lb(isLower));
            dxu        = tau*(ub(isUpper) - x(isUpper));
            x(isLower) = lb(isLower) + dxl;
            x(isUpper) = ub(isUpper) - dxu;
        end


    end

end