classdef DualUpdater_NullSpace < DualUpdater

    properties (Access = private)
       NSmeritFunc
       constraint
       cost
       problem
       options
    end

    methods (Access = public)

        function obj = DualUpdater_NullSpace(cParams)
            obj.init(cParams)
            obj.constraint   = cParams.constraint;
            obj.cost         = cParams.cost;
            %obj.NSmeritFunc  = cParams.NSmeritFunc;
        end

        function updateDualVariable(obj)
            obj.computeDualProblemParameters();
            obj.computeDualProblemOptions();
            PROBLEM         = obj.problem;
            PROBLEM.options = obj.options;
            mu = quadprog(PROBLEM);
            obj.dualVariable.value = mu;
        end

    end

    methods (Access = private)

        function computeDualProblemParameters(obj)
            A = obj.constraint.gradient;
            b = obj.cost.gradient;
            prob.H      = A'*A;
            prob.f      = b'*A;
            prob.A      = [];
            prob.b      = [];
            prob.Aeq    = [];
            prob.beq    = [];
            prob.lb     = -inf;
            prob.ub     = inf;
            prob.x0     = zeros(length(prob.H),1);
            prob.solver = 'quadprog';
            obj.problem = prob;
        end

        function computeDualProblemOptions(obj)
            opts = optimoptions("quadprog");
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
