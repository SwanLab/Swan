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
            obj.dualVariable.restart();
            lambda = obj.dualVariable.value;            
            fref = obj.computeFeasibleDesignVariable(lambda);
            
            isLB = false;
            isUB = false;
            i = -15;
            pow = 1.1;
           
            while ~isLB && ~isUB && i < 1000
                lLB = lambda - pow^(i);
                fLB = obj.computeFeasibleDesignVariable(lLB);
                
                
                lUB = lambda + pow^(i);
                fUB = obj.computeFeasibleDesignVariable(lUB);
                
                isLB = fLB*fref < 0;
                isUB = fUB*fref < 0;
                i = i + 1;
            end
            
             if isLB
                 lUB = lambda + pow^(i-2);
             else
                 lLB = lambda - pow^(i-2);
             end
             
             if i > 49
                a = 0; 
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
            obj.constraint.computeFunctionAndGradient();
            fval = obj.constraint.value;
        end
        
        function updateDesignVariable(obj)
            obj.unconstrainedOptimizer.hasConverged = false;
            obj.unconstrainedOptimizer.compute();
            obj.unconstrainedOptimizer.updateConvergenceParams();
        end
        
    end
    
end