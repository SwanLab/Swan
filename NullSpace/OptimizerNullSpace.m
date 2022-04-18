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
        tol = 1e-4
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
%                 obj.checkConvergence();
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
            obj.costOld    = obj.cost.value;
            obj.designVariable.updateOld();
            obj.lambda     = 0;
        end

        function obj = update(obj)
            obj.calculateInitialStep();
            obj.acceptableStep   = false;
            obj.lineSearchTrials = 0;
            x0   = obj.designVariable.value;
            obj.saveOldValues(x0);
            mOld = obj.computeMeritFunction(x0);            
            while ~obj.acceptableStep
                obj.updateDualDirect();
                x = obj.updatePrimal();
                obj.checkStep(x,x0,mOld);
            end
            obj.updateOldValues(x);
        end

        function obj = calculateInitialStep(obj)
            obj.cost.computeFunctionAndGradient();
            x       = obj.designVariable.value;
            g       = obj.cost.gradient;
            if isempty(obj.tau)
                obj.tau = 0.1*sqrt(norm(x)/norm(g));
            else
                obj.tau = 3*obj.tau;
            end
        end

        function obj = updateDualDirect(obj)
            obj.constraint.computeFunctionAndGradient();
            obj.cost.computeFunctionAndGradient();
            x = obj.designVariable.value;
            A = obj.constraint.gradient;
            c = x-obj.tau*obj.cost.gradient;
            b = A'*x - obj.constraint.value ;
            l = (A'*A)\((c'*A)' - b);
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
            lb     = obj.lowerBound;
            ub     = obj.upperBound;
            t      = obj.tau;
            A      = obj.constraint.gradient;
            b      = obj.cost.gradient;
            l      = obj.lambda;
            x      = obj.designVariable.value;            
            c      = x - t*b;
%             S      = (A'*A)^-1;
%             nu     = -t*S*A'*b - l*A'*S*A;
%             xN     = c - A*nu;
            xN     = c - A*l;
%             aJ = 1;
%             aC = 0.1;
%             I  = ones(1,size(A,2));
%             eJ = (I - A*S*A')*b;
%             eC = S*A;
%             xN = x - t*(aJ*eJ + aC*l*eC);
            x  = min(ub,max(xN,lb));
%             alphaC = 0.01*obj.tau;            
%             xN     = x0 - t*b - alphaC*A*(A'*A)^(+1)*l;
        end

        function checkStep(obj,x,x0,mOld)
            mNew = obj.computeMeritFunction(x);
            if mNew < mOld
                obj.acceptableStep = true;
            elseif obj.tau < 1e-10
                obj.acceptableStep = true;
                obj.tau = 0.1;
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
            x0        = obj.oldDesignVariable;
            l         = obj.lambda;
            C         = obj.constraint.value;
%             DC        = obj.constraint.gradient';
            J         = obj.cost.value;
            DJ        = obj.cost.gradient';
            t         = obj.tau;
%             alphaJ    = 1;
%             alphaC    = 1;
%             S         = (DC*DC')^-1;
%             mF        = alphaJ*(J + l*C) + alphaC/2*C'*S*C;
            mF        = J + l'*C; %+ (DJ + l*C)*(x-x0) + 1/(2*t)*norm(x-x0)^2;
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
            obj.nIter = obj.nIter + 1;
        end

        function obj = checkConvergence(obj)
           if abs(obj.oldCost - obj.cost.value) < obj.tol && max(obj.constraint.value) < 0
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