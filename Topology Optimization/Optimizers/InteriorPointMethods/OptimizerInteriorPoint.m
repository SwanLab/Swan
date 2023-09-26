classdef OptimizerInteriorPoint < Optimizer

    properties (Access = private)
        slack
        lowerX
        upperX
        constraintLowerBound
        constraintUpperBound
        tol = 1e-5
        lowerSlack
        upperSlack
        lowerZ
        upperZ
        barrierTau
        lineSearch
        lineSearchTrials
        error
        H
        hessian
        dx, ds 
        dlam
        dzL, dzU
        zLa, zUa
        xNew, sNew
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
            obj.dualUpdater.compute(obj.lowerZ,obj.upperZ);
        end

        function moveVariablesToFeasible(obj)
            s.x      = obj.designVariable.value;
            s.upperX = obj.upperX;
            s.lowerX = obj.lowerX;
            x        = obj.replaceOutOfBoundsVariables(s);
            obj.designVariable.update(x);
            s.x       = obj.slack;
            s.upperX  = obj.upperSlack;
            s.lowerX  = obj.lowerSlack;
            obj.slack = obj.replaceOutOfBoundsVariables(s);
        end

        function computeDualVariableBounds(obj)
            mu = obj.baseVariables.mu;
            x  = obj.designVariable.value';
            lx = obj.lowerX;
            ux = obj.upperX;
            s  = obj.slack';
            ls = obj.lowerSlack;
            us = obj.upperSlack;
            zL = mu./(x-lx);
            zL = [zL,mu./(s-ls)];
            zU = mu./(ux-x);
            zU = [zU,mu./(us-s)];
            obj.lowerZ = zL;
            obj.upperZ = zU;
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
            s.lowerBound   = obj.lowerX';
            s.upperBound   = obj.upperX';
            obj.cost.value = c - mu*obj.computeLogPenaltyTerm(s);
        end

        function penalizeDueToSlackVariable(obj)
            c              = obj.cost.value;
            mu             = obj.baseVariables.mu;
            s.field        = obj.slack;
            s.lowerBound   = obj.constraintLowerBound(obj.isInequality());
            s.upperBound   = obj.constraintUpperBound(obj.isInequality());
            obj.cost.value = c - mu*obj.computeLogPenaltyTerm(s);
        end

        function itIs = isInequality(obj)
            itIs = obj.constraintUpperBound > obj.constraintLowerBound;
        end

        function updateCostGradient(obj)
            nS = size(obj.slack,2);
            obj.cost.gradient = [obj.cost.gradient' zeros(1,nS)];
        end

        function updateConstraint(obj)
            g                      = obj.constraint.value;
            gLB                    = obj.constraintLowerBound(:,1);
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
            x   = obj.designVariable.value;
            g   = -obj.dx;
            x   = obj.primalUpdater.update(g,x);
        end

        function updatePrimalVariables(obj)
            alphaPrimal = obj.primalUpdater.tau;
            if obj.nSlack >= 1
                obj.sNew      = obj.slack + alphaPrimal * obj.ds';
            else
                obj.sNew      = [];
            end
        end

        function updateDual(obj)
            alphaDual = obj.primalUpdater.tau;
            obj.dualUpdater.update(alphaDual,obj.dlam);
            obj.zLa   = obj.lowerZ + alphaDual*obj.dzL';
            obj.zUa   = obj.upperZ + alphaDual*obj.dzU';           
        end

        function pushVariables(obj)
            s.field      = obj.xNew;
            s.lowerBound = obj.lowerX;
            s.upperBound = obj.upperX;
            s.tau        = obj.barrierTau;
            obj.xNew     = obj.pushVarsFunction(s);
            s.field      = obj.slack;
            s.lowerBound = obj.lowerSlack;
            s.upperBound = obj.upperSlack;
            obj.sNew     = obj.pushVarsFunction(s);
            obj.zLa(obj.zLa<0) = obj.barrierTau*(obj.lowerZ(obj.zLa<0));
            obj.zUa(obj.zUa<0) = obj.barrierTau*(obj.upperZ(obj.zUa<0));
        end

        function updateWithAcceptance(obj)
            obj.designVariable.update(obj.xNew);
            obj.slack  = obj.sNew;
            obj.lowerZ = obj.zLa;
            obj.upperZ = obj.zUa;
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
            s.lowerBounds.X  = obj.lowerX;
            s.lowerBounds.Z  = obj.lowerZ;
            s.lowerBounds.s  = obj.lowerSlack;
            s.upperBounds.X  = obj.upperX;
            s.upperBounds.Z  = obj.upperZ;
            s.upperBounds.s  = obj.upperSlack;
            dirs             = IPMDirectionComputer(s);
            dirs.compute();
            obj.dx           = dirs.gradients.dx;
            obj.ds           = dirs.gradients.ds;
            obj.dlam         = dirs.gradients.dlam;
            obj.dzL          = dirs.gradients.dzL;
            obj.dzU          = dirs.gradients.dzU;
            obj.H            = dirs.updatedHessian;
        end
        
        function init(obj,cParams)
            obj.upperX                 = cParams.uncOptimizerSettings.ub;
            obj.lowerX                 = cParams.uncOptimizerSettings.lb;
            obj.cost                   = cParams.cost;
            obj.constraint             = cParams.constraint;
            obj.designVariable         = cParams.designVar;
            obj.dualVariable           = cParams.dualVariable;
            obj.incrementalScheme      = cParams.incrementalScheme;
            obj.maxIter                = cParams.maxIter;
            obj.nIter                  = 0;
            obj.constraintCase         = cParams.constraintCase;
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
                    obj.constraintLowerBound(i,:) = -inf;
                    obj.constraintUpperBound(i,:) = 0;
                else
                    obj.constraintLowerBound(i,:) = 0;
                    obj.constraintUpperBound(i,:) = 0;
                end
            end
        end

        function correctUpperAndLowerBounds(obj)
            nnode      = obj.designVariable.mesh.nnodes;
            obj.lowerX = (obj.lowerX-1e-6)*ones(1,nnode);
            obj.upperX = (obj.upperX+1e-6)*ones(1,nnode);
        end

        function computeSlackVariables(obj)
            obj.lowerSlack = obj.constraintLowerBound(obj.isInequality());
            obj.upperSlack = obj.constraintUpperBound(obj.isInequality());
            obj.slack      = zeros(size(obj.lowerSlack));
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
            sMax      = 100;
            nConstr   = obj.constraint.nSF;
            nnode     = obj.designVariable.mesh.nnodes;
            e         = ones(nnode + obj.nSlack,1);
            sD        = max(sMax,(sum(abs(obj.dualVariable.value))+sum(abs(obj.lowerZ))+sum(abs(obj.upperZ)))/(nConstr+2*(nnode+obj.nSlack)));
            sC        = max(sMax,(sum(abs(obj.upperZ))+sum(abs(obj.lowerZ)))/(2*(nnode+obj.nSlack)));
            part(1)   = max(abs(obj.cost.gradient' + obj.constraint.gradient*obj.dualVariable.value' - obj.lowerZ' + obj.upperZ'))/sD;
            part(2)   = max(abs(obj.constraint.value));
            part(3)   = max(abs(diag([obj.designVariable.value'-obj.lowerX obj.slack-obj.lowerSlack])*diag(obj.lowerZ)*e))/sC;
            part(4)   = max(abs(diag([obj.upperX-obj.designVariable.value' obj.upperSlack-obj.slack])*diag(obj.upperZ)*e))/sC;
            obj.error = max(part);
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

        function x = replaceOutOfBoundsVariables(s)
            x          = s.x;
            upperX     = s.upperX;
            lowerX     = s.lowerX;
            isLower    = x<=lowerX;
            isUpper    = x>=upperX;
            x(isLower) = min(upperX(isLower),lowerX(isLower)+1e-2);
            x(isUpper) = max(lowerX(isUpper),upperX(isUpper)-1e-2);
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