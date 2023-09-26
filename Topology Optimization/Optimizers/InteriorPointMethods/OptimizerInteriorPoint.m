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
        e
        nConstr
        alphaPrimal
        alphaDual
        barrierTau

        lineSearch
        lineSearchTrials
        e_mu
        LHS
        RHS
        explicitSol
        invDiagdLX, invDiagdUX
        lowerSigma, upperSigma
        H
        hessian
        dx, ds 
        dlam, lam
        dzL, dzU
        zLa, zUa
        xa, sa
        hasConverged
        hasFinished
        incrementalScheme
        nSlack
        nX
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
            obj.computeVariableConstraintMultipliers();
            obj.updateCC();
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

        function computeVariableConstraintMultipliers(obj)
            obj.barrierTau = min(obj.baseVariables.tau_max,100*obj.baseVariables.mu);
            zL = obj.baseVariables.mu./(obj.designVariable.value' - obj.lowerX);
            zL = [zL,obj.baseVariables.mu./(obj.slack' - obj.lowerSlack)];
            zU = obj.baseVariables.mu./(obj.upperX - obj.designVariable.value');
            zU = [zU,obj.baseVariables.mu./(obj.upperSlack - obj.slack')];

            obj.lowerZ = zL;
            obj.upperZ = zU;
        end

        function updateCost(obj)
             x    = obj.designVariable.value;
             lb   = obj.lowerX';
             ub   = obj.upperX';
             mu   = obj.baseVariables.mu;
             c    = obj.cost.value;             
             cNew = c - mu*sum(log(x-lb) + log(ub-x));

             sl  = obj.slack;
             gLB = obj.constraintLowerBound(obj.isInequality());
             gUB = obj.constraintUpperBound(obj.isInequality());
             cNew      = cNew - sum(mu*(log(sl - gLB) + log(gUB - sl)));

             obj.cost.value = cNew;
        end

        function itIs = isInequality(obj)
             itIs = obj.constraintUpperBound > obj.constraintLowerBound;
        end

        function updateCostGradient(obj)
            nS = size(obj.slack,2);
            obj.cost.gradient = [obj.cost.gradient' zeros(1,nS)];
        end

        function updateConstraint(obj)
            g   = obj.constraint.value;
            gUB = obj.constraintUpperBound(:,1);
            gLB = obj.constraintLowerBound(:,1);
            sl = obj.slack;
            areBoundsEqual = gUB == gLB;
            g(areBoundsEqual)    = g(areBoundsEqual) - gLB(areBoundsEqual);
            g(~areBoundsEqual)   = g(~areBoundsEqual) - sl;
            obj.constraint.value = g;
        end

        function updateConstraintGradient(obj)
            mC = size(obj.constraint.gradient,1);
            nC = size(obj.constraint.gradient,2);
            isUpper = obj.constraintUpperBound > obj.constraintLowerBound;
            dk = length(find(isUpper==true));
            i  = mC+1:mC+dk;
            obj.constraint.gradient(i,nC) = -1;
        end

        function computeHessian(obj)
            hess      = obj.hessian(1:obj.nX,1:obj.nX); 
            deltaX    = obj.designVariable.value - obj.oldDesignVariable;
            deltaCost = obj.cost.gradient - obj.oldCostGradient;
            if deltaCost == zeros(size(deltaCost))
                obj.hessian = zeros(obj.nX,obj.nX);
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
            J  = obj.cost.value;
            g  = obj.constraint.value;
            meritF = J + obj.baseVariables.nu*sum(abs(g));
            obj.updateCostGradient();
            obj.updateConstraint();
            obj.updateConstraintGradient();
        end

        function updateCC(obj)
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
                obj.xa = obj.updatePrimal();
                obj.updatePrimalVariables();
                obj.updateDual();
                obj.pushVariables();
                obj.checkStep(obj.xa,x0);
            end
        end

        function obj = calculateInitialStep(obj)
            x          = obj.designVariable.value;
            DJ         = obj.cost.gradient;
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
            obj.alphaPrimal = obj.primalUpdater.tau;
            if (obj.nSlack >= 1)
                obj.sa      = obj.slack + obj.alphaPrimal * obj.ds';
            else
                obj.sa      = [];
            end
        end

        function updateDual(obj)
            obj.alphaDual = obj.primalUpdater.tau;
            obj.dualUpdater.update(obj.alphaDual,obj.dlam);
            obj.zLa       = obj.lowerZ + obj.alphaDual * obj.dzL';
            obj.zUa       = obj.upperZ + obj.alphaDual * obj.dzU';           
        end

        function pushVariables(obj)
            isLower = obj.xa<obj.lowerX;
            isUpper = obj.xa>obj.upperX;
            dxl     = obj.barrierTau*(obj.designVariable.value(isLower) - obj.lowerX(isLower));
            dxu     = obj.barrierTau*(obj.upperX(isUpper) - obj.designVariable.value(isUpper));
            obj.xa(isLower) = obj.lowerX(isLower) + dxl;
            obj.xa(isUpper) = obj.upperX(isUpper) - dxu;
            obj.zLa(obj.zLa<0) = obj.barrierTau*(obj.lowerZ(obj.zLa<0));
            obj.zUa(obj.zUa<0) = obj.barrierTau*(obj.upperZ(obj.zUa<0));
            isLower = obj.sa<obj.lowerSlack;
            isUpper = obj.sa>obj.upperSlack;
            dsl     = obj.barrierTau*(obj.slack(isLower) - obj.lowerSlack(isLower));
            dsu     = obj.barrierTau*(obj.upperSlack(isUpper) - obj.slack(isUpper));
            obj.sa(isLower) = obj.lowerSlack(isLower) + dsl;
            obj.sa(isUpper) = obj.upperSlack(isUpper) - dsu;
        end

        function updateWithAcceptance(obj)
            obj.designVariable.update(obj.xa);
            obj.slack  = obj.sa;
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
            costPrediction(1)         = obj.cost.gradient*[obj.dx;obj.ds];
            costPrediction(2)         = max(0,0.5*[obj.dx;obj.ds]'*obj.H*[obj.dx;obj.ds]);
            costDecreasePrediction    = sum(costPrediction);
            theta                     = sum(abs(obj.constraint.value));
            rho                       = 0.1;
            nuUpdated                 = costDecreasePrediction/((1 - rho) * theta);
            obj.baseVariables.nu      = max(1,min(1000,nuUpdated));
        end

        function computeOptimizerDirections(obj)
            s.cost           = obj.cost;
            s.constraint     = obj.constraint;
            s.designVariable = obj.designVariable;
            s.slack          = obj.slack;
            s.dualVariable   = obj.dualVariable;
            s.baseVariables  = obj.baseVariables;
            s.hessian        = obj.hessian;
            s.nConstr        = obj.nConstr;
            s.nX             = obj.nX;
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
            obj.H            = dirs.H;
        end
        
        function init(obj,cParams)
            obj.upperX                 = cParams.uncOptimizerSettings.ub;
            obj.lowerX                 = cParams.uncOptimizerSettings.lb;
            obj.cost                   = cParams.cost;
            obj.constraint             = cParams.constraint;
            obj.designVariable         = cParams.designVar;
            obj.dualVariable           = cParams.dualVariable;
            obj.incrementalScheme      = cParams.incrementalScheme;
            obj.nX                     = length(obj.designVariable.value);
            obj.maxIter                = cParams.maxIter;
            obj.nIter                  = 0;
            obj.constraintCase         = cParams.constraintCase;
            obj.cost.computeFunctionAndGradient();
            obj.hessian = eye(obj.nX);
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
            obj.computeUpperAndLowerBounds();
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

        function computeUpperAndLowerBounds(obj)
            n = obj.nX;
            obj.lowerX = (-1e-6)*ones(1,n);
            obj.upperX = (1+1e-6)*ones(1,n);
        end

        function computeSlackVariables(obj)
            obj.nConstr   = max(size(obj.constraintLowerBound));
            isUpper = obj.constraintUpperBound > obj.constraintLowerBound;
            slPos   = 1:length(find(isUpper==true));
            obj.slack(slPos)      = 0;
            obj.lowerSlack(slPos) = obj.constraintLowerBound(isUpper);
            obj.upperSlack(slPos) = obj.constraintUpperBound(isUpper);
            obj.nSlack            = max(size(obj.slack));
            obj.updateConstraint();
            obj.nConstr = max(size(obj.constraint.value));
            obj.e = ones(obj.nX + obj.nSlack,1);
            if (obj.baseVariables.slack_init)
                obj.slack(slPos) = obj.constraint.value(isUpper);
            else
                obj.slack(slPos) = 0.01;
            end
            obj.alphaPrimal = 1.0;
            obj.alphaDual   = 1.0;
        end

        function checkConvergence(obj)
            obj.cost.computeFunctionAndGradient();
            obj.computeHessian();
            obj.constraint.computeFunctionAndGradient();
            obj.updateCC();
            obj.computeError();
            if obj.e_mu <= obj.tol
                obj.hasConverged = true;
            end
        end
        function computeError(obj)
            s_max    = 100; % > 1
            s_d      = max(s_max,(sum(abs(obj.dualVariable.value))+sum(abs(obj.lowerZ))+sum(abs(obj.upperZ)))/(obj.nConstr+2*(obj.nX+obj.nSlack)));
            s_c      = max(s_max,(sum(abs(obj.upperZ))+sum(abs(obj.lowerZ)))/(2*(obj.nX+obj.nSlack)));
            part(1)  = max(abs(obj.cost.gradient' + obj.constraint.gradient*obj.dualVariable.value' - obj.lowerZ' + obj.upperZ'))/s_d;
            part(2)  = max(abs(obj.constraint.value));
            part(3)  = max(abs(diag([obj.designVariable.value'-obj.lowerX obj.slack-obj.lowerSlack])*diag(obj.lowerZ)*obj.e))/s_c;
            part(4)  = max(abs(diag([obj.upperX-obj.designVariable.value' obj.upperSlack-obj.slack])*diag(obj.upperZ)*obj.e))/s_c;
            obj.e_mu = max(part);
        end

        function checkNewBarrierProblem(obj)
            k_mu                     = obj.baseVariables.k_mu;
            if (obj.e_mu < k_mu * obj.baseVariables.mu)
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

    end

end