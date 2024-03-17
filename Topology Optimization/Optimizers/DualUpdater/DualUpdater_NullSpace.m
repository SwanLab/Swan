classdef DualUpdater_NullSpace < handle

    properties (Access = private)
       constraint
       cost
       designVariable
       dualVariable
       options
       constraintCase
       nConstr
       dualOld
       isInequality
       position
    end

    methods (Access = public)

        function obj = DualUpdater_NullSpace(cParams)
            obj.init(cParams);
            obj.defineConstraintCases();
            obj.computeDualProblemOptions();
        end

        function defineConstraintCases(obj)
            k = 1;
            for i = 1:obj.nConstr
                switch obj.constraintCase{i}
                    case {'EQUALITY'}
                       obj.isInequality  = false;
                    case {'INEQUALITY'}
                        obj.isInequality = true;
                        obj.position(k)  = i;
                        k                = k + 1;
                end
            end
        end

        function prob = computeDualBounds(obj)
            tol     = inf;
            prob.lb = -tol*ones(obj.nConstr,1);
            prob.ub = tol*ones(obj.nConstr,1);
            if obj.isInequality
                p          = obj.position;
                prob.lb(p) = 0;
            end
        end

        function update(obj,eta,ub,lb)
            s.prob = obj.computeDualBounds();
            s.eta  = eta;
            s.ub   = ub;
            s.lb   = lb;
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
            obj.designVariable = cParams.designVariable;
            obj.constraintCase = cParams.constraintCase;
            obj.dualVariable   = cParams.dualVariable;
            obj.dualOld        = obj.dualVariable.value;
            obj.nConstr        = length(cParams.dualVariable.value);
        end

        function computeQuadraticProblem(obj,s)
%             lJ = obj.computeDualNullSpace();
%             lG = obj.computeDualRangeSpace(s);
            l  = obj.computeDualGlobal(s);
            obj.dualVariable.value = l;
        end

        function lJ = computeDualNullSpace(obj)
            Dg              = obj.constraint.gradient;
            DJ              = obj.cost.gradient;
            problem.H       = Dg'*Dg;
            problem.f       = Dg'*DJ;
            problem.solver  = 'quadprog';
            problem.options = obj.options;
            lJ              = quadprog(problem);
        end

        function lG = computeDualRangeSpace(obj,s)
            eta             = s.eta;
            Dg              = obj.constraint.gradient;
            g               = obj.constraint.value;
            problem.H       = Dg'*Dg;
            problem.f       = -eta*g;
            problem.solver  = 'quadprog';
            problem.options = obj.options;
            lG              = quadprog(problem);
        end

        function l = computeDualGlobal(obj,s)
            % Main point in next meeting 18/03/2024:
            eta     = s.eta;
            problem = s.prob;
            ub      = s.ub;
            lb      = s.lb;
            g       = obj.constraint.value;
            Dg      = obj.constraint.gradient;
            DJ      = obj.cost.gradient;
            l       = obj.dualVariable.value;
            if sum(l)~=0
                dMref = DJ+Dg*l;
                x = obj.designVariable.fun.fValues;
                switch class(obj.designVariable)
                    case 'LevelSet'
                        isUBActive = find(x<=0 & dMref<0);
                        isLBActive = find(x>0 & dMref>0);
                    otherwise
                        isUBActive = find(abs(x-ub)<=1e-8 & dMref<0);
                        isLBActive = find(abs(x-lb)<=1e-8 & dMref>0);
                end
                noAct = [isUBActive;isLBActive];
            else
                noAct = [];
            end
            Dg(noAct)       = 0;
            DJ(noAct)       = 0;
            problem.H       = Dg'*Dg;
            problem.f       = Dg'*DJ-eta*g;
            problem.solver  = 'quadprog';
            problem.options = obj.options;
            l               = quadprog(problem);
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

end