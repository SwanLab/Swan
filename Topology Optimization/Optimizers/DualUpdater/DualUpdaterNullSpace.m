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
       primalUpdater
       acceptableStep
    end

    methods (Access = public)
        function obj = DualUpdaterNullSpace(cParams)
            obj.init(cParams);
            obj.defineConstraintCases();
            obj.computeDualProblemOptions();
            obj.createPrimalUpdater();
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
                t     = primalUpdater.boxConstraints.refTau;
                s.lUB = primalUpdater.boxConstraints.lUB/t;
                s.lLB = primalUpdater.boxConstraints.lLB/t;
            end
            obj.computeQuadraticProblem(s);
        end

        function reset(obj)
            obj.dualVariable.value = obj.dualOld;
        end

        function updateOld(obj)
            obj.dualOld = obj.dualVariable.fun.fValues;
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.constraintCase = cParams.constraintCase;
            obj.dualVariable   = cParams.dualVariable;
            obj.dualOld        = obj.dualVariable.fun.fValues;
            obj.nConstr        = length(cParams.dualVariable.fun.fValues);
        end

        function is = isEquality(obj,i)
            switch obj.constraintCase{i}
                case 'EQUALITY'
                    is = true;
                case 'INEQUALITY'
                    is = false;
            end
        end

        function is = isInequality(obj,i)
            switch obj.constraintCase{i}
                case 'EQUALITY'
                    is = false;
                case 'INEQUALITY'
                    is = true;
            end        
        end

        function computeQuadraticProblem(obj,s)
            eta      = s.eta;
            lUB      = s.lUB;
            lLB      = s.lLB;
            problem  = s.prob;
            g        = obj.constraint.value;
            Dg       = obj.constraint.gradient;
            isActive = obj.checkComplementaryKKT(g);
            problem.lb(~isActive) = [];
            problem.ub(~isActive) = [];
            g  = g(isActive);
            Dg = Dg(:,isActive);
            DJ = obj.cost.gradient;
            l  = zeros(obj.nConstr,1);
            if ~isempty(g)
                problem.H       = Dg'*Dg;
                problem.f       = Dg'*(DJ+lUB-lLB)-eta*g;
                problem.solver  = 'quadprog';
                problem.options = obj.options;
                l(isActive)     = quadprog(problem);
            end
            obj.dualVariable.fun.fValues = l;

            %obj.solveWithProjectedGradient(s);
        end

        function isActive = checkComplementaryKKT(obj,g)
            isActive = true(obj.nConstr,1);
            for i = 1:obj.nConstr
                if obj.isInequality(i) && g(i)<-1e-2
                    isActive(i) = false;
                end
            end
        end

        function solveWithProjectedGradient(obj,s)
            eta   = s.eta;
            lUB   = s.lUB;
            lLB   = s.lLB;
            l     = obj.dualVariable;
            g     = obj.constraint.value;
            Dg    = obj.constraint.gradient;
            DJ    = obj.cost.gradient;
            A     = Dg'*Dg;
            b     = eta*g-Dg'*(DJ+lUB-lLB);
            Df    = A*l.fun.fValues-b;
            obj.primalUpdater.tau = Df'*Df/(Df'*A*Df);
            delta = inf;
            while delta>1e-8
                lValk              = l.fun.fValues;
                r                  = b-A*lValk;
                obj.acceptableStep = false;
                while (~obj.acceptableStep)
                    l     = obj.primalUpdater.update(-r,l);
                    lValk1   = l.fun.fValues;
                    obj.checkStep(A,b,lValk,lValk1);
                end
                delta = norm(lValk1-lValk);
                obj.primalUpdater.increaseStepLength(1.2);
            end
            obj.dualVariable = l;
        end

        function checkStep(obj,A,b,lOld,l)
            Mold = obj.computeMeritFunction(A,b,lOld);
            M    = obj.computeMeritFunction(A,b,l);
            if M<Mold+1e-6
                obj.acceptableStep = true;
            else
                obj.primalUpdater.decreaseStepLength();
                obj.dualVariable.update(lOld);
            end
        end
        
        function M = computeMeritFunction(obj,A,b,l)
            lLB = obj.primalUpdater.boxConstraints.lLB;
            q   = 0.5*l'*A*l-l'*b;
            M   = q-sum(lLB);
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

        function createPrimalUpdater(obj)
            s                 = obj.computeDualBounds();
            p                 = ProjectedGradient(s);
            obj.primalUpdater = p;
        end
    end
end