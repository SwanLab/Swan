classdef ConstraintProjector < handle

    properties (Access = private)
        problem
        cost
        constraint
        dualVariable
        designVariable
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
            x0      = obj.designVariable.fun.fValues;
            tolCons = 1e-5;
            obj.tau = obj.primalUpdater.tau;
            tolCons = 1e-2*tolCons;
            lambda  = obj.dualVariable.fun.fValues;
            fref    = obj.computeFeasibleDesignVariable(x0,lambda);
            if abs(fref) > tolCons
                obj.computeBounds(x0);
                obj.problem.x0 = [obj.lambdaLB obj.lambdaUB];
                obj.problem.options = optimset(obj.problem.options,'TolX',tolCons);
                obj.problem.objective = @(lambda) obj.computeFeasibleDesignVariable(x0,lambda);
                fzero(obj.problem);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams,s)
            obj.cost             = cParams.cost;
            obj.constraint       = cParams.constraint;
            obj.designVariable   = cParams.designVariable;
            obj.dualVariable     = cParams.dualVariable;
            obj.primalUpdater    = s.primalUpdater;
        end

        function defineProblem(obj)
            obj.problem.solver    = 'fzero';
            obj.problem.options   = optimset(@fzero);
        end

        function computeBounds(obj,x0)
            lambda = obj.dualVariable.fun.fValues;
            fref   = obj.computeFeasibleDesignVariable(x0,lambda);
            isLB   = false;
            isUB   = false;
            i      = -15;
            pow    = 3;
            while ~isLB && ~isUB && i < 1000
                lLB  = lambda - pow^(i);
                fLB  = obj.computeFeasibleDesignVariable(x0,lLB);
                lUB  = lambda + pow^(i);
                fUB  = obj.computeFeasibleDesignVariable(x0,lUB);
                isLB = fLB*fref < 0;
                isUB = fUB*fref < 0;
                i    = i + 1;
            end
            if isLB
                 lUB = lambda + pow^(i-2);
            else
                 lLB = lambda - pow^(i-2);
            end
            obj.lambdaLB = lLB;
            obj.lambdaUB = lUB;
        end

        function fval = computeFeasibleDesignVariable(obj,x0,lambda)
            obj.designVariable.update(x0);
            d = obj.designVariable;
            obj.dualVariable.fun.fValues = lambda;
            obj.updatePrimal();
            obj.constraint.computeFunctionAndGradient(d);
            fval = obj.constraint.value;
        end

        function x = updatePrimal(obj)
            Dg = obj.constraint.gradient;
            DJ = obj.cost.gradient;
            l  = obj.dualVariable.fun.fValues;
            x  = obj.designVariable;
            g  = DJ + l*Dg;
            x  = obj.primalUpdater.update(g,x);
            obj.designVariable = x;
        end

    end

end
