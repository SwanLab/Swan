classdef InteriorPointMethodsSolver < Optimizer
    properties (GetAccess = public, SetAccess = protected)
        type = 'IPM';
    end
    properties (Access = private)
        slack
        lowerX
        upperX
        lowerLinearBound
        upperLinearBound
        tol = 1e-8
        lowerSlack
        upperSlack
        lowerZ
        upperZ
        e
        m
        alphaPrimal
        alphaDual
        barrierTau

        tau
        lineSearch
        lineSearchTrials
        isposdef
        e_mu
        lowerDiagonalZ
        upperDiagonalZ
        diagonaldL
        diagonaldU
        dL
        dU
        LHS
        RHS
        explicitSol
        invDiagdL, invDiagdU
        lowerSigma, upperSigma
        H
        hessian
        alphaDualMax
        alphaPrimalMax
        dx, ds 
        dlam, lam
        dzL, dzU
        zLa, zUa
        xa, sa
        accept
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
        function obj = InteriorPointMethodsSolver(cParams)
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
            %obj.finalSolTest();
        end
    end

    methods(Access = private)
        function previousComputations(obj)
            obj.moveVariablesToFeasible();
            obj.computeVariableConstraintMultiplyers();
            obj.updateCC();
            obj.dualUpdater.compute(obj.lowerZ,obj.upperZ);
        end

        function moveVariablesToFeasible(obj)
            % move x into feasible region and off boundary initially
            for i = 1:obj.nX
                if (obj.designVariable.value(i) <= obj.lowerX(i) || obj.designVariable.value(i) >= obj.upperX(i))
                    fprintf(1,'Moving x0 to interior region\n')
                    break
                end
            end
            % Move x variables to be feasible
            for i = 1:obj.nX
                if (obj.designVariable.value(i) <= obj.lowerX(i))
                    obj.designVariable.value(i) = min(obj.upperX(i),obj.lowerX(i)+1e-2);
                end
                if (obj.designVariable.value(i) >= obj.upperX(i))
                    obj.designVariable.value(i) = max(obj.lowerX(i),obj.upperX(i)-1e-2);
                end
            end
            % Move slack variables to be feasible
            for i = 1:obj.nSlack
                if (obj.slack(i) <= obj.lowerSlack(i))
                    obj.slack(i) = min(obj.upperSlack(i),obj.lowerSlack(i)+1e-2);
                end
                if (obj.slack(i) >= obj.upperSlack(i))
                    obj.slack(i) = max(obj.lowerSlack(i),obj.upperSlack(i)-1e-2);
                end
            end
        end

        function computeVariableConstraintMultiplyers(obj)
            obj.barrierTau = min(obj.baseVariables.tau_max,100*obj.baseVariables.mu);
            % variable constraint multipliers
            % zL*(x-xL) = mu  =>  zL = mu / (x-xL)
            % zL(1:n) = 1;
            for i = 1:obj.nX
                zL(i) = obj.baseVariables.mu / (obj.designVariable.value(i) - obj.lowerX(i));
            end
            for i = 1:obj.nSlack
                zL(obj.nX + i) = obj.baseVariables.mu / (obj.slack(i) - obj.lowerSlack(i));
            end
            % zU*(xU-x) = mu  =>  zU = mu / (xU-x)
            for i = 1:obj.nX
                zU(i) = obj.baseVariables.mu / (obj.upperX(i) - obj.designVariable.value(i));
            end
            for i = 1:obj.nSlack
                zU(obj.nX+i) = obj.baseVariables.mu / (obj.upperSlack(i) - obj.slack(i));
            end
            obj.lowerZ = zL;
            obj.upperZ = zU;
        end

        function updateCost(obj)
             for i = 1:obj.nX              
                obj.cost.value = obj.cost.value - obj.baseVariables.mu * (log(obj.designVariable.value(i) - obj.lowerX(i)) + log(obj.upperX(i) - obj.designVariable.value(i)));
             end
             j = 0;
             for i = 1:obj.m
                if(obj.upperLinearBound(i) > obj.lowerLinearBound(i))
                   j = j + 1;
                   obj.cost.value = obj.cost.value - obj.baseVariables.mu * (log(obj.slack(j) - obj.lowerLinearBound(i)) + log(obj.upperLinearBound(i) - obj.slack(j)));
                end
             end
        end

        function updateCostGradient(obj)
            nS = size(obj.slack,2);
            obj.cost.gradient = [obj.cost.gradient' zeros(1,nS)];
        end

        function updateConstraint(obj)
            j = 0;
            for i = 1:size(obj.constraint.value,1) %changed 2 for 1
                if (obj.upperLinearBound(i) == obj.lowerLinearBound(i))
                    % equality constant
                    obj.constraint.value(i) = obj.constraint.value(i) - obj.lowerLinearBound(i);
                else
                    % inequality slack
                    j = j + 1;
                    obj.constraint.value(i) = obj.constraint.value(i) - obj.slack(j);
                end
            end
        end

        function updateConstraintGradient(obj)
            mC = size(obj.constraint.gradient,1);
            nC = size(obj.constraint.gradient,2);
            k = 0;
            for i = 1:mC
                if(obj.upperLinearBound(i) > obj.lowerLinearBound(i))
                    k = k + 1;
                    obj.constraint.gradient(i,nC+k) = -1;
                end
            end
        end

        function hessianComputer(obj)
            hess = obj.hessian(1:obj.nX,1:obj.nX); 
            deltaX = obj.designVariable.value - obj.oldDesignVariable;
            deltaCost = obj.cost.gradient - obj.oldCostGradient;
            if deltaCost == zeros(size(deltaCost))
                obj.hessian = zeros(obj.nX,obj.nX);
            else
                costFrac = (deltaCost*deltaCost')/(deltaCost'*deltaX);
                xFrac = (hess*deltaX*(hess*deltaX)')/(deltaX'*hess*deltaX);
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
            %obj.updateCC();
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
            obj.sigmaComputer();
            obj.updateHessian();            % -> updateHessianSize()
            obj.linearSystemSolver();
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
            x       = obj.designVariable.value;
            if obj.nIter == 0
                factor = 1;
                obj.primalUpdater.computeFirstStepLength(x,x,factor);
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
                obj.sa = obj.slack + obj.alphaPrimal * obj.ds';
            else
                obj.sa = [];
            end
        end

        function updateDual(obj)
            obj.alphaDual = obj.primalUpdater.tau;
            obj.lam = obj.dualVariable.value;
            obj.dualUpdater.update(obj.alphaDual,obj.dlam);
            obj.zLa = obj.lowerZ + obj.alphaDual * obj.dzL';
            obj.zUa = obj.upperZ + obj.alphaDual * obj.dzU';
        end

        function pushVariables(obj)
            % clipping (this should already be arranged by alpha_max)
            % push away from the boundary with tau
            for i = 1:obj.nX
                if(obj.xa(i) < obj.lowerX(i))
                    obj.xa(i) = obj.lowerX(i) + obj.barrierTau*(obj.designVariable.value(i) - obj.lowerX(i));
                end
                if(obj.xa(i) > obj.upperX(i))
                    obj.xa(i) = obj.upperX(i) - obj.barrierTau*(obj.upperX(i) - obj.designVariable.value(i));
                end
                if(obj.zLa(i) < 0)
                    obj.zLa(i) = obj.barrierTau*(obj.lowerZ(i));
                end
                if(obj.zUa(i) < 0)
                    obj.zUa(i) = obj.barrierTau*(obj.upperZ(i));
                end
            end
            for i = 1:obj.nSlack
                if(obj.sa(i) < obj.lowerSlack(i))
                    obj.sa(i) = obj.lowerSlack(i) + obj.barrierTau*(obj.slack(i) - obj.lowerSlack(i));

                end
                if(obj.sa(i) > obj.upperSlack(i))
                    obj.sa(i) = obj.upperSlack(i) - obj.barrierTau*(obj.upperSlack(i) - obj.slack(i));
                end
            end
        end

        function updateWithAcceptance(obj)
            obj.designVariable.update(obj.xa);
            obj.slack = obj.sa;
            obj.lowerZ = obj.zLa;
            obj.upperZ = obj.zUa;
        end

        function updateHessian(obj)
            for i = 1:obj.nSlack
                obj.hessian(obj.nX+i,obj.nX+i) = 0.0;
            end
            obj.H = obj.hessian + obj.lowerSigma + obj.upperSigma;
        end

        function obj = saveOldValues(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.oldCost            = obj.cost.value;
            obj.oldCostGradient = obj.cost.gradient;
            obj.designVariable.updateOld;
            obj.oldDesignVariable = obj.designVariable.value;
        end

        function checkStep(obj,x,x0)
            mNew = obj.computeMeritFunction(x);
            pred_phi_1 = obj.cost.gradient*[obj.dx;obj.ds];
            pred_phi_2 = max(0,0.5*[obj.dx;obj.ds]'*obj.H*[obj.dx;obj.ds]);
            pred_phi_decrease = pred_phi_1 + pred_phi_2;
            theta = sum(abs(obj.constraint.value));
            rho = 0.1;
            nuUpdated = pred_phi_decrease / ((1 - rho) * theta);
            obj.baseVariables.nu = max(1,min(1000,nuUpdated));
            predictUpdated = -obj.alphaPrimal*obj.cost.gradient*[obj.dx;obj.ds] - 0.5*obj.alphaPrimal^2*[obj.dx;obj.ds]'*obj.H*[obj.dx;obj.ds] + obj.baseVariables.nu*(norm(obj.constraint.value',1) - norm(obj.constraint.value' + obj.alphaPrimal*obj.constraint.gradient'*[obj.dx;obj.ds],1));
            reduced = obj.mOld - mNew;
            eta = 0.2;
            if  reduced >= eta*predictUpdated
%             if mNew < obj.mOld
                obj.acceptableStep = true;
%                 obj.primalUpdater.tau = 1;
%                 obj.alphaPrimal = 1;
%                 obj.alphaDual = 1;
                obj.updateWithAcceptance();
                obj.meritNew = mNew;
            elseif obj.primalUpdater.isTooSmall()
                warning('Convergence could not be achieved (step length too small)')
                obj.acceptableStep = true;
                obj.alphaPrimal = 1;
                obj.alphaDual = 1;
                obj.meritNew = mNew; % Provisional value
            else
                obj.primalUpdater.decreaseStepLength();
                obj.designVariable.update(x0);
                obj.lineSearchTrials = obj.lineSearchTrials + 1;
            end
        end

        function sigmaComputer(obj)
            obj.lowerDiagonalZ = diag(obj.lowerZ);
            obj.upperDiagonalZ = diag(obj.upperZ);
            % sigmas
            obj.dL = [obj.designVariable.value'-obj.lowerX obj.slack-obj.lowerSlack];
            obj.dU = [obj.upperX-obj.designVariable.value' obj.upperSlack-obj.slack];
            obj.diagonaldL = diag(obj.dL);
            obj.diagonaldU = diag(obj.dU);
            invdML = zeros(obj.nX+obj.nSlack,obj.nX+obj.nSlack);
            invdMU = zeros(obj.nX+obj.nSlack,obj.nX+obj.nSlack);
            sigL = zeros(obj.nX+obj.nSlack,obj.nX+obj.nSlack);
            sigU = zeros(obj.nX+obj.nSlack,obj.nX+obj.nSlack);
            for i = 1:obj.nX+obj.nSlack
                invdML(i,i) = 1 / obj.diagonaldL(i,i); 
                invdMU(i,i) = 1 / obj.diagonaldU(i,i); 
                sigL(i,i) = obj.lowerDiagonalZ(i,i) / obj.diagonaldL(i,i);
                sigU(i,i) = obj.upperDiagonalZ(i,i) / obj.diagonaldU(i,i);
             end
             obj.invDiagdL = invdML;
             obj.invDiagdU = invdMU;
             obj.lowerSigma = sigL;
             obj.upperSigma = sigU;
        end

        function linearSystemSolver(obj)
            % compute search direction, solving Ax = b
            obj.computeLHS();
            obj.computeRHS();
            obj.explicitSol = -pinv(full(obj.LHS))*obj.RHS;
            obj.searchZDirection();
        end

        function searchZDirection(obj)
            obj.dx = obj.explicitSol(1:obj.nX,1);
            obj.ds = obj.explicitSol(obj.nX+1:obj.nX+obj.nSlack,1);
            obj.dlam = obj.explicitSol(obj.nX + obj.nSlack + 1:obj.nX + obj.nSlack + obj.m,1);
            obj.dzL = obj.baseVariables.mu*obj.invDiagdL*obj.e - obj.lowerZ' - obj.lowerSigma*[obj.dx; obj.ds];
            obj.dzU = obj.baseVariables.mu*obj.invDiagdU*obj.e - obj.upperZ' + obj.upperSigma*[obj.dx; obj.ds];
        end

        function computeLHS(obj)
            s.H = obj.H;
            s.m = obj.m;
            s.constraint = obj.constraint;
            s.funcType = 'Symmetric';
            cLHS = LHSComputer.create(s);
            cLHS.compute();
            obj.LHS = cLHS.LHS;
        end

        function computeRHS(obj)
            s.baseVariables = obj.baseVariables;
            s.e = obj.e;
            s.nX = obj.nX;
            s.nSlack = obj.nSlack;
            s.m = obj.m;
            s.cost = obj.cost;
            s.constraint = obj.constraint;
            s.lambda = obj.dualVariable.value;
            s.lowerZ = obj.lowerZ;
            s.upperZ = obj.upperZ;
            s.invDiagdL = obj.invDiagdL;
            s.invDiagdU = obj.invDiagdU;
            s.funcType = 'Symmetric';
            cRHS = RHSComputer.create(s);
            cRHS.compute();
            obj.RHS = cRHS.RHS;
        end
        
        function init(obj,cParams)
            obj.upperX = cParams.uncOptimizerSettings.ub;
            obj.lowerX = cParams.uncOptimizerSettings.lb;
            obj.cost                   = cParams.cost;
            obj.constraint             = cParams.constraint;
            obj.designVariable         = cParams.designVar;
            obj.dualVariable           = cParams.dualVariable;
            obj.incrementalScheme      = cParams.incrementalScheme;
            obj.nX                     = length(obj.designVariable.value);
            obj.maxIter                = cParams.maxIter;
            obj.hasConverged           = false;
            obj.nIter                  = 0;
            obj.constraintCase         = cParams.constraintCase;
            obj.cost.computeFunctionAndGradient();
            obj.hessian = eye(obj.nX);
            obj.constraint.computeFunctionAndGradient();
            obj.loadIPMVariables();
        end

        function loadIPMVariables(obj)
            obj.baseVariables.nu = 0.1;
            obj.baseVariables.mu = 10;
            obj.baseVariables.slack_init = true;
            % update tau
            obj.baseVariables.tau_max = 0.01;
            obj.baseVariables.k_mu  = 0.2;
            obj.computeLinearBounds();
            obj.computeSlackVariables();
            obj.computeUpperAndLowerBounds();
            obj.baseVariables.idebug = 0;
        end

        function computeLinearBounds(obj)
            for i = 1:length(obj.constraintCase)
                if strcmp(obj.constraintCase{i},'INEQUALITY')
                    for j = 1:obj.nX
                        obj.lowerLinearBound(j) = -inf;
                        obj.upperLinearBound(j) = 0;
                    end
                else
                    for j = 1:obj.nX
                        obj.lowerLinearBound(j) = 0;
                        obj.upperLinearBound(j) = 0;
                    end
                end
            end
        end

        function computeUpperAndLowerBounds(obj)
            for i = 1:obj.nX
                obj.lowerX(i) = -inf;
                obj.upperX(i) = inf;
            end
        end

        function computeSlackVariables(obj)
            % add slack for inequality constraints only
            k = 0;
            obj.m = max(size(obj.lowerLinearBound));
            obj.slack= [];
            obj.lowerSlack = [];
            obj.upperSlack = [];
            for i = 1:obj.m
                if(obj.upperLinearBound(i) > obj.lowerLinearBound(i))
                    k = k + 1;
                    obj.slack(k) = 0;
                    obj.lowerSlack(k) = obj.lowerLinearBound(i);
                    obj.upperSlack(k) = obj.upperLinearBound(i);
                end
            end
            obj.nSlack = max(size(obj.slack));
            obj.m = max(size(obj.constraint.value));
            obj.e = ones(obj.nX + obj.nSlack,1);
            k = 0;
            for i = 1:obj.m
                if(obj.upperLinearBound(i) > obj.lowerLinearBound(i))
                    k = k + 1;
                    if (obj.baseVariables.slack_init)
                    obj.slack(k) = obj.constraint.value(i);
                    else
                    obj.slack(k) = 0.01;
                    end
                end
            end
            obj.alphaPrimal = 1.0;
            obj.alphaDual = 1.0;
        end

        function checkConvergence(obj)
            x = obj.designVariable.value;
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.hessianComputer();
            obj.constraint.computeFunctionAndGradient();
            obj.updateCC();
            s_max = 100; % > 1
            s_d = max(s_max,(sum(abs(obj.dualVariable.value))+sum(abs(obj.lowerZ))+sum(abs(obj.upperZ)))/(obj.m+2*(obj.nX+obj.nSlack)));
            s_c = max(s_max,(sum(abs(obj.upperZ))+sum(abs(obj.lowerZ)))/(2*(obj.nX+obj.nSlack)));
            part(1) = max(abs(obj.cost.gradient' + obj.constraint.gradient*obj.dualVariable.value' - obj.lowerZ' + obj.upperZ'))/s_d;
            part(2) = max(abs(obj.constraint.value));
            part(3) = max(abs(diag([obj.designVariable.value'-obj.lowerX obj.slack-obj.lowerSlack])*diag(obj.lowerZ)*obj.e))/s_c;
            part(4) = max(abs(diag([obj.upperX-obj.designVariable.value' obj.upperSlack-obj.slack])*diag(obj.upperZ)*obj.e))/s_c;
            obj.e_mu = max(part);
            if obj.e_mu <= obj.tol
                %obj.finalSolTest();
                obj.hasConverged = true;
            end
        end

        function checkNewBarrierProblem(obj)
            k_mu = obj.baseVariables.k_mu;
            if (obj.e_mu < k_mu * obj.baseVariables.mu)
                th_mu = 1.5;
                obj.baseVariables.mu = max(obj.tol/10,min(k_mu*obj.baseVariables.mu,obj.baseVariables.mu^th_mu));
                obj.barrierTau = min(obj.baseVariables.tau_max,100*obj.baseVariables.mu);
            end
        end

        function obj = updateMonitoring(obj)
            s.nIter            = obj.nIter;
            s.tau              = obj.primalUpdater.tau;
            obj.lineSearch = 'PROJECTED GRADIENT';
            s.lineSearch       = obj.lineSearch;
            s.lineSearchTrials = obj.lineSearchTrials;
            s.oldCost          = obj.oldCost;
            s.hasFinished      = obj.hasFinished;
            s.meritNew         = obj.meritNew;
            obj.outputFunction.monitoring.compute(s);
            if obj.hasConverged == false
                obj.printIteration();
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
            iStep = obj.incrementalScheme.iStep;
            nStep = obj.incrementalScheme.nSteps;
            itHas = obj.nIter >= obj.maxIter*(iStep/nStep);
        end

        function printIteration(obj)
            theta = sum(abs(obj.constraint.value));
            du = sum(abs(obj.cost.gradient' + obj.constraint.gradient'*obj.dualVariable.value' - obj.lowerZ' + obj.upperZ'));
            logmu = log10(obj.baseVariables.mu);
            if (mod(obj.nIter,15)==1)
            fprintf(1,'   Iter       Merit   Objective   log10(mu)        Pcvg        Dcvg    alpha_pr    alpha_du\n');
            end
            fprintf(1,'  %5i %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n',obj.nIter,obj.meritNew,obj.cost.value,logmu,theta,du,obj.alphaPrimal,obj.alphaDual);
        end
    end
end