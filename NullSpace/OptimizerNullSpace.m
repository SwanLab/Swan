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
        lambda
        incrementalScheme
        hasFinished
        mOld
    end

    methods (Access = public) 
        
        function obj = OptimizerNullSpace(cParams)
            obj.initOptimizer(cParams);
            obj.init(cParams);
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
            obj.lambda = 0;
        end

        function obj = update(obj)
            if obj.nIter == 10
                a = 0;
            end
            x0   = obj.designVariable.value;
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
            l       = obj.lambda;
            DJ      = obj.cost.gradient;
            Dg      = obj.constraint.gradient;
            aJ      = 1;
            DmF     = aJ*(DJ + l*Dg);
%             obj.tau = 1*sqrt(norm(DmF)/norm(x));
            if obj.nIter == 0
                obj.tau = 1*sqrt(norm(DmF)/norm(x));
            else
                obj.tau = 1.05*obj.tau;
            end
            obj.tau = min(obj.tau,1.2);
        end

        function obj = updateDualDirect(obj)
            obj.constraint.computeFunctionAndGradient();
            obj.cost.computeFunctionAndGradient();
            DJ = obj.cost.gradient;
            Dg = obj.constraint.gradient;
            g  = obj.constraint.value;
            S  = (Dg'*Dg)^-1;
            aC = 1;
            aJ = 1;
            t  = obj.tau;
            l  = aC/aJ*S*(g - t*Dg'*DJ);
            obj.lambda = l;
        end

%         function obj = updateDualQuadProg(obj)
%             obj.constraint.computeFunctionAndGradient();
%             obj.cost.computeFunctionAndGradient();
%             obj.computeDualProblemParameters();
%             obj.computeDualProblemOptions();
%             PROBLEM         = obj.problem;
%             PROBLEM.options = obj.options;
%             l = quadprog(PROBLEM);
%             obj.lambda = l;
%         end
% 
%          function computeDualProblemParameters(obj)
%             A = obj.constraint.gradient;
%             b = obj.cost.gradient;
%             c = obj.constraint.value;
%             prob.H      = A'*A;
%             prob.f      = c - obj.tau*b'*A;
%             prob.A      = [];
%             prob.b      = [];
%             prob.Aeq    = [];
%             prob.beq    = [];
%             prob.lb     = -inf;
%             prob.ub     = inf;
%             prob.x0     = zeros(length(prob.H),1);
%             prob.solver = 'quadprog';
%             obj.problem = prob;
%         end
% 
%         function computeDualProblemOptions(obj)
%             opts = optimoptions("quadprog");
%             opts = struct( ...
%                 'Algorithm','interior-point-convex', ...
%                 'Diagnostics','off', ...
%                 'Display','none', ...
%                 'HessMult',[], ...
%                 'MaxIter',1e3, ...
%                 'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
%                 'PrecondBandWidth',0, ...
%                 'ProblemdefOptions', struct, ...
%                 'TolCon',1e-5, ...
%                 'TolFun',[], ...
%                 'TolFunValue', [], ...
%                 'TolPCG',0.1, ...
%                 'TolX',100*eps, ...
%                 'TypicalX','ones(numberOfVariables,1)', ...
%                 'LinearSolver', 'auto', ...
%                 'ObjectiveLimit', -1e20 ...
%                 );
%             obj.options = opts;
%         end
% 
        function x = updatePrimal(obj)
            lb      = obj.lowerBound;
            ub      = obj.upperBound;
            t       = obj.tau;
            Dg      = obj.constraint.gradient;
            g       = obj.constraint.value;
            DJ      = obj.cost.gradient;
            l       = obj.lambda;
            x       = obj.designVariable.value;
            S       = (Dg'*Dg)^-1;
            aJ  = 1;
            aC  = 1;
            dAJ     = aJ*(DJ + l*Dg);
            dAC     = aC*Dg'*S*Dg;
            dx      = -t*(dAJ);
            xN      = x + dx;
            x       = min(ub,max(xN,lb));
        end

        function checkStep(obj,x,x0)
            mNew = obj.computeMeritFunction(x);

            if obj.nIter == 0 && mNew == obj.mOld
                obj.acceptableStep = true;
                obj.updateDualDirect();
            end
            if mNew < obj.mOld
                obj.acceptableStep = true;
                obj.updateDualDirect();
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
            DJ     = obj.cost.gradient;
            J      = obj.cost.value;
            Dg     = obj.constraint.gradient;
            g      = obj.constraint.value;
            l      = obj.lambda;
            S      = (Dg'*Dg)^-1;
            aJ     = 1;
            aC     = 0.1;
            AJ     = aJ*(J + l*g);
            AC     = aC/2*g'*S*g;
            mF     = AJ;
        end

        function obj = saveOldValues(obj,x)
            obj.designVariable.update(x);
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.oldCost            = obj.cost.value;
            obj.oldDesignVariable  = x;
            obj.dualVariable.value = obj.lambda;
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