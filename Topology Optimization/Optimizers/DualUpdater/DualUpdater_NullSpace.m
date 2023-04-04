classdef DualUpdater_NullSpace < handle

    properties (Access = private)
       constraint
       cost
       designVariable
       dualVariable
       problem
       options
       constraintCase
       nConstr
       dualOld
       isInequality
       position
       parameter
    end

    methods (Access = public)

        function obj = DualUpdater_NullSpace(cParams)
            obj.init(cParams);
            obj.defineConstraintCases();
            obj.defineRangeStepParameterValue(cParams);
        end

        function defineConstraintCases(obj)
            obj.isInequality = false;
            k                = 1;
            for i = 1:obj.nConstr
                switch obj.constraintCase{i}
                    case {'EQUALITY'}
                       
                    case {'INEQUALITY'}
                        obj.isInequality = true;
                        obj.position(k)  = i;
                        k                = k + 1;
                end
            end
        end

        function defineRangeStepParameterValue(obj,cParams)
            switch cParams.optimizerNames.primal
                case {'PROJECTED GRADIENT','HAMILTON JACOBI'}
                    obj.parameter = inf;
                case {'SLERP'}
                    obj.parameter = 5;
                otherwise

            end
        end

        function update(obj)
            if obj.isInequality
                obj.computeQuadraticProblem();
            else
                obj.computeDirectDual();
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

        function computeDirectDual(obj)
            obj.constraint.computeFunctionAndGradient();
            obj.cost.computeFunctionAndGradient();
            DJ = obj.cost.gradient;
            Dh = obj.constraint.gradient;
            h  = obj.constraint.value;
            S  = (Dh'*Dh)^-1;
            aJ = 1;
            aC = 1;
            f  = obj.parameter;
            f2 = 1;
            AC = min(f,f2*aC/aJ*S*h);
            AJ = -aC/aJ*S*Dh'*DJ;
            l  = AC + AJ;
            obj.dualVariable.value = l;
        end

        function computeQuadraticProblem(obj)
            obj.constraint.computeFunctionAndGradient();
            obj.cost.computeFunctionAndGradient();
            obj.computeDualProblemParameters();
            obj.computeDualProblemOptions();
            PROBLEM         = obj.problem;
            PROBLEM.options = obj.options;
            l = quadprog(PROBLEM);
            obj.dualVariable.value = l;
        end

         function computeDualProblemParameters(obj)
            Dg          = obj.constraint.gradient;
            DJ          = obj.cost.gradient;
            g           = obj.constraint.value;
            f           = obj.parameter;
            c           = min(f*(Dg'*Dg),g);
            i           = length(g);
            p           = obj.position;
            prob.H      = Dg'*Dg;
            prob.f      = -c + Dg'*DJ;
            prob.A      = [];
            prob.b      = [];
            prob.Aeq    = [];
            prob.beq    = [];
            prob.lb     = -inf*ones(i,1);
            prob.lb(p)  = 0;
            prob.ub     = inf*ones(length(g),1);
            prob.x0     = zeros(length(prob.H),1);
            prob.solver = 'quadprog';
            obj.problem = prob;
        end

        function computeDualProblemOptions(obj)
            %opts = optimoptions('quadprog');
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
