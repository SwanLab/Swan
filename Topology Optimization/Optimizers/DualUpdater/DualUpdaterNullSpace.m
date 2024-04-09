classdef DualUpdaterNullSpace < handle

    properties (Access = private)
       constraint
       cost
       dualVariable
       options
       constraintCase
       nConstr
       dualOld
       position
    end

    methods (Access = public)
        function obj = DualUpdaterNullSpace(cParams)
            obj.init(cParams);
            obj.defineConstraintCases();
            obj.computeDualProblemOptions();
        end

        function defineConstraintCases(obj)
            k = 1;
            for i = 1:obj.nConstr
                switch obj.constraintCase{i}
                    case {'INEQUALITY'}
                        obj.position(k) = i;
                        k               = k + 1;
                end
            end
        end

        function prob = computeDualBounds(obj)
            tol        = inf;
            prob.lb    = -tol*ones(obj.nConstr,1);
            prob.ub    = tol*ones(obj.nConstr,1);
            p          = obj.position;
            prob.lb(p) = 0;
        end

        function update(obj,eta,primalUpdater)
            s.prob = obj.computeDualBounds();
            s.eta  = eta;
            if isempty(primalUpdater.tau)
                s.lUB = 0;
                s.lLB = 0;
            else
                s.lUB = primalUpdater.lUB/primalUpdater.tau;
                s.lLB = primalUpdater.lLB/primalUpdater.tau;
            end
            obj.computeQuadraticProblem(s);
        end

        function reset(obj)
            obj.dualVariable.value = obj.dualOld;
        end

        function updateOld(obj)
            obj.dualOld = obj.dualVariable.value;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.constraintCase = cParams.constraintCase;
            obj.dualVariable   = cParams.dualVariable;
            obj.dualOld        = obj.dualVariable.value;
            obj.nConstr        = length(cParams.dualVariable.value);
        end

        function computeQuadraticProblem(obj,s)
            eta             = s.eta;
            lUB             = s.lUB;
            lLB             = s.lLB;
            problem         = s.prob;
            g               = obj.constraint.value;
            Dg              = obj.constraint.gradient;
            DJ              = obj.cost.gradient;
            problem.H       = Dg'*Dg;
            problem.f       = Dg'*(DJ+lUB-lLB)-eta*g;
            problem.solver  = 'quadprog';
            problem.options = obj.options;
            l               = quadprog(problem);
            obj.dualVariable.value = l;
        end

        function computeDualProblemOptions(obj)
            opts = struct( ...
                'Algorithm','interior-point-convex', ...
                'Diagnostics','off', ...
                'Display','none', ...
                'HessMult',[], ...
                'MaxIter',1e3, ...
                'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
                'PrecondBandWidth',0, ...
                'ProblemdefOptions', struct, ...
                'TolCon',1e-5, ...
                'TolFun',[], ...
                'TolFunValue', [], ...
                'TolPCG',0.1, ...
                'TolX',100*eps, ...
                'TypicalX','ones(numberOfVariables,1)', ...
                'LinearSolver', 'auto', ...
                'ObjectiveLimit', -1e20 ...
                );
            obj.options = opts;
        end
    end
end