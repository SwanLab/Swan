classdef DualUpdater_NullSpace < handle

    properties (Access = private)
       NSmeritFunc
       constraint
       cost
       designVariable
       dualVariable
       problem
       options
       constraintCase
       tau
       nConstr
       constrTol
       dualOld
    end

    methods (Access = public)

        function obj = DualUpdater_NullSpace(cParams)
            obj.init(cParams);
        end

        function update(obj)
            for i = 1:obj.nConstr
                switch obj.constraintCase{i}
                    case {'EQUALITY'}
                        obj.computeDirectDual(i);
                    case {'INEQUALITY'}
                        obj.computeQuadraticProblem(i);

                end
            end
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
            obj.designVariable = cParams.designVar;
            obj.constraintCase = cParams.constraintCase;
            obj.dualVariable   = cParams.dualVariable;
            obj.dualOld        = obj.dualVariable.value;
            obj.nConstr        = cParams.constraint.nSF;
        end

        function computeDirectDual(obj,i)
            obj.constraint.computeFunctionAndGradient();
            obj.cost.computeFunctionAndGradient();
            DJ = obj.cost.gradient;
            Dh = obj.constraint.gradient(:,i);
            h  = obj.constraint.value(i);
            S  = (Dh'*Dh)^-1;
            aJ = 1;
            aC = 1;
            l  = aC/aJ*S*(h - 1*Dh'*DJ);
            obj.dualVariable.value(i) = l;
        end

        function computeQuadraticProblem(obj,i)
            obj.constraint.computeFunctionAndGradient();
            obj.cost.computeFunctionAndGradient();
            obj.computeDualProblemParameters(i);
            obj.computeDualProblemOptions();
            PROBLEM         = obj.problem;
            PROBLEM.options = obj.options;
            l = quadprog(PROBLEM);
            obj.dualVariable.value(i) = l;
        end

         function computeDualProblemParameters(obj,i)
            Dg = obj.constraint.gradient(:,i);
            DJ = obj.cost.gradient;
            g  = obj.constraint.value(i);
            t  = 1;%obj.tau;
            prob.H      = Dg'*Dg;
            prob.f      = -g + t*Dg'*DJ;
            prob.A      = [];
            prob.b      = [];
            prob.Aeq    = [];
            prob.beq    = [];
            prob.lb     = 0;
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
