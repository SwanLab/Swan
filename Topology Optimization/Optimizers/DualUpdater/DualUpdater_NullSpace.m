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

    properties (Access = public)
        t
        aGMax
        aG

        lGtrialPl
        lGmaxPl
        lGPl
        lJPl
        lPl
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
            S  = (Dh'*Dh)^-1;
            g  = obj.constraint.value;
            lJ = -S*Dh'*DJ;
            lGtrial = obj.aG*S*g/obj.t;
            lGmax   = obj.aGMax*max(abs(DJ+lJ*Dh))/max(abs(Dh));
            lG = obj.projectLambdaG(lGtrial,lGmax);
            l  = lG + lJ;
            obj.dualVariable.value = l;

            obj.lGtrialPl = lGtrial;
            obj.lGmaxPl = lGmax;
            obj.lGPl = lG;
            obj.lJPl = lJ;
            obj.lPl = l;
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
            c           = min(f*(Dg'*Dg),g);% min(f*(Dg'*DJ),g)
            i           = length(g);
            p           = obj.position;
            prob.H      = Dg'*Dg;
            prob.f      = -c + Dg'*DJ;
            prob.A      = [];
            prob.b      = [];
            prob.Aeq    = [];
            prob.beq    = [];

            tol = 15; % inf

            prob.lb     = -tol*ones(i,1);
            prob.lb(p)  = 0;
            prob.ub     = tol*ones(length(g),1);
            prob.x0     = zeros(length(prob.H),1);
            prob.solver = 'quadprog';
            obj.problem = prob;
        end

        function computeDualProblemOptions(obj)
%             opts = optimoptions("quadprog");
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

    methods (Static,Access = private)

        function lG = projectLambdaG(lG,Delta)
            lG = max(-Delta,min(Delta,lG));
        end

    end


end
