classdef ConstraintProjector < handle

    properties (Access = private)
        problem
        cost
        constraint
        dualVariable
        designVariable
        targetParameters
        primalUpdater
        lambdaUB
        lambdaLB
        tau
    end

    methods (Access = public)

        function obj = ConstraintProjector(cParams,s)
            obj.init(cParams,s);
            obj.defineProblem();
        end

        function project(obj)
            obj.tau = obj.primalUpdater.tau;
            tolCons = 1e-2*obj.targetParameters.constr_tol;
            lambda  = obj.dualVariable.value;
            fref    = obj.computeFeasibleDesignVariable(lambda);
            if abs(fref) > tolCons
                obj.computeBounds();
                obj.problem.x0 = [obj.lambdaLB obj.lambdaUB];
                obj.problem.options = optimset(obj.problem.options,'TolX',tolCons);
                fzero(obj.problem);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams,s)
            obj.cost             = cParams.cost;
            obj.constraint       = cParams.constraint;
            obj.designVariable   = cParams.designVar;
            obj.dualVariable     = cParams.dualVariable;
            obj.targetParameters = cParams.targetParameters;
            obj.primalUpdater    = s.primalUpdater;
        end

        function defineProblem(obj)
            obj.problem.solver    = 'fzero';
            obj.problem.options   = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.computeFeasibleDesignVariable(lambda);
        end

        function computeBounds(obj)
            lambda = obj.dualVariable.value;
            fref   = obj.computeFeasibleDesignVariable(lambda);
            isLB   = false;
            isUB   = false;
            i      = -15;
            pow    = 1.1; %1.1
            while ~isLB && ~isUB && i < 1000 % i < 1000
                lLB  = lambda - pow^(i);
                fLB  = obj.computeFeasibleDesignVariable(lLB);
                lUB  = lambda + pow^(i);
                fUB  = obj.computeFeasibleDesignVariable(lUB);
                isLB = fLB*fref < 0;
                isUB = fUB*fref < 0;
                i    = i + 1;
            end
            display(i)
            if isLB
                 lUB = lambda + pow^(i-2);
            else
                 lLB = lambda - pow^(i-2);
            end
            obj.lambdaLB = lLB;
            obj.lambdaUB = lUB;
        end

        function fval = computeFeasibleDesignVariable(obj,lambda)
            obj.designVariable.restart();
            obj.dualVariable.value = lambda;
            obj.updatePrimal();
            obj.constraint.computeFunctionAndGradient();
            fval = obj.constraint.value;
        end

        function x = updatePrimal(obj)
            Dg = obj.constraint.gradient;
            DJ = obj.cost.gradient;
            l  = obj.dualVariable.value;
            x  = obj.designVariable.value;
            g  = DJ + l*Dg;
            x  = obj.primalUpdater.update(g,x);
            obj.designVariable.update(x);
        end

    end

end
