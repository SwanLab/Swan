classdef ConstraintProjector < handle
    
    properties (Access = private)
        problem
        lagrangian
        cost
        constraint
        dualVariable
        designVariable
        targetParameters
        unconstrainedOptimizer
        lambdaUB
        lambdaLB
    end
    
    methods (Access = public)
        
        function obj = ConstraintProjector(cParams)
            obj.cost       = cParams.cost;
            obj.constraint = cParams.constraint;
            obj.designVariable = cParams.designVariable;
            obj.dualVariable = cParams.dualVariable;
            obj.targetParameters = cParams.targetParameters;
            obj.lagrangian = cParams.lagrangian;
            obj.unconstrainedOptimizer = cParams.unconstrainedOptimizer;
            obj.defineProblem();
        end
        
        function project(obj)
            tolCons = 1e-2*obj.targetParameters.constr_tol;  
            lambda = obj.dualVariable.value;
            fref = obj.computeFeasibleDesignVariable(lambda);
            if abs(fref) > tolCons              
                obj.computeBounds();    
                obj.problem.x0 = [obj.lambdaLB obj.lambdaUB];
                obj.problem.options = optimset(obj.problem.options,'TolX',tolCons);
                fzero(obj.problem);
            end
        end
        
    end
    
    methods (Access = private)
        
        function defineProblem(obj)
            obj.problem.solver = 'fzero';
            obj.problem.options = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.computeFeasibleDesignVariable(lambda);
        end
        
        function computeBounds(obj)
            lambda = obj.dualVariable.value;            
            fref = obj.computeFeasibleDesignVariable(lambda);
            

            lLB = lambda - 10^-4;
            fLB = obj.computeFeasibleDesignVariable(lLB);
            
            if fLB < fref
                if fref < 0
                    isUB2change = true;
                else
                    isUB2change = false;
                end
            else
                if fref < 0
                    isUB2change = false;
                else
                    isUB2change = true;
                end
            end
            
            fB = fref;
            i = -4;
            while fref*fB > 0
                if isUB2change
                    lambdaB = lambda + 10^i;
                else
                    lambdaB = lambda - 10^i;
                end
                 fB = obj.computeFeasibleDesignVariable(lambdaB);

                i = i + 1;
            end
            
            if isUB2change
                lLB = lambda;
                lUB = lambdaB;
            else
                lLB = lambdaB;
                lUB = lambda;
            end
            obj.lambdaLB = lLB;
            obj.lambdaUB = lUB;
            
        end
        
        function fval = computeFeasibleDesignVariable(obj,lambda)
            obj.designVariable.restart();
            obj.constraint.restart();
            obj.dualVariable.value = lambda;
            obj.lagrangian.computeGradient();
            obj.updateDesignVariable();
            obj.constraint.computeCostAndGradient();
            fval = obj.constraint.value;
        end
        
        function updateDesignVariable(obj)
            obj.unconstrainedOptimizer.hasConverged = false;
            obj.unconstrainedOptimizer.compute();
            obj.unconstrainedOptimizer.updateConvergenceParams();
        end
        
    end
    
end