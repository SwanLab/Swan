classdef DualUpdaterNullSpace < handle

    properties (Access = private)
        cost
        constraint
        constraintCase
        nConstr
    end

    properties (Access = private)
        position
        l0
    end

    methods (Access = public)
        function obj = DualUpdaterNullSpace(cParams)
            obj.init(cParams);
            obj.defineConstraintCases();
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

        function l = update(obj,eta,lUB,lLB)
            s.prob = obj.computeDualBounds();
            s.eta  = eta;
            s.lUB = lUB;
            s.lLB = lLB;
            l = obj.computeQuadraticProblem(s);
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.cost           = cParams.cost;
            obj.constraint     = cParams.constraint;
            obj.constraintCase = cParams.constraintCase;
            obj.nConstr        = length(cParams.constraintCase);
        end

        function is = isInequality(obj,i)
            switch obj.constraintCase{i}
                case 'EQUALITY'
                    is = false;
                case 'INEQUALITY'
                    is = true;
            end
        end

        function l = computeQuadraticProblem(obj,s)
            eta           = s.eta;
            xBoxUB        = s.lUB;
            xBoxLB        = s.lLB;
            lb            = s.prob.lb;
            ub            = s.prob.ub;
            g             = obj.constraint.value;
            Dg            = obj.constraint.gradient;
            isActive      = obj.checkComplementaryKKT(g);
            lb(~isActive) = [];
            ub(~isActive) = [];
            g  = g(isActive);
            Dg = Dg(:,isActive);
            DJ = obj.cost.gradient;
            l  = zeros(obj.nConstr,1);
            if ~isempty(g)
                H           = Dg'*Dg;
                f           = Dg'*(DJ+xBoxUB-xBoxLB)-eta*g;
                l(isActive) = obj.solve(H,f,lb,ub,isActive);
            end
            obj.l0 = l;
        end

        function isActive = checkComplementaryKKT(obj,g)
            isActive = true(obj.nConstr,1);
            for i = 1:obj.nConstr
                if obj.isInequality(i) && g(i)<-1e-2
                    isActive(i) = false;
                end
            end
        end

        function l = solve(obj,H,f,lb,ub,isActive)
            s.cost           = obj.createCost(H,f);
            s.designVariable = obj.createDesignVariable(isActive);
            s.monitoring     = false;
            s.lb             = lb;
            s.ub             = ub;
            s.maxIter        = 1000;
            opt              = OptimizerProjectedGradient(s);
            opt.solveProblem();
            l = s.designVariable.fun.fValues;
        end

        function d = createDesignVariable(obj,isActive)
            g    = obj.constraint.value;
            s.x0 = zeros(size(g(isActive)));
            d    = DesignVariableAcademic(s);
        end
    end

    methods (Static, Access = private)
        function c = createCost(H,f)
            sA.cF = @(x) 0.5*x'*H*x + f'*x;
            sA.gF = @(x) H*x + f;
            s.shapeFunctions{1} = AcademicCost(sA);
            s.weights           = 1;
            s.Msmooth           = 1;
            c                   = Cost(s);
        end
    end
end