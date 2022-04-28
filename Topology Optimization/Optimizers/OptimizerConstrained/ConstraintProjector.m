classdef ConstraintProjector < handle
    
    properties (Access = private)
        problem
        cost
        constraint
        dualVariable
        designVariable
        targetParameters
%         unconstrainedOptimizer
        lambdaUB
        lambdaLB
        x
        tau
        lowerBound
        upperBound
    end
    
    methods (Access = public)
        
        function obj = ConstraintProjector(cParams)
            init(cParams);
            obj.defineProblem();
        end
        
        function project(obj,x0)
            obj.x   = x0;
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

        function init(obj,cParams)
            obj.cost             = cParams.cost;
            obj.constraint       = cParams.constraint;
            obj.designVariable   = cParams.designVariable;
            obj.dualVariable     = cParams.dualVariable;
            obj.targetParameters = cParams.targetParameters;
            obj.tau              = cParams.tau;
            obj.upperBound       = cParams.upperBound;
            obj.lowerBound       = cParams.lowerBound;
        end
        
        function defineProblem(obj)
            obj.problem.solver    = 'fzero';
            obj.problem.options   = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.computeFeasibleDesignVariable(lambda);
        end
        
        function computeBounds(obj)
            obj.dualVariable.restart();
            lambda = obj.dualVariable.value;            
            fref   = obj.computeFeasibleDesignVariable(lambda);            
            isLB   = false;
            isUB   = false;
            i      = -15;
            pow    = 1.1;           
            while ~isLB && ~isUB && i < 1000
                lLB  = lambda - pow^(i);
                fLB  = obj.computeFeasibleDesignVariable(lLB);                               
                lUB  = lambda + pow^(i);
                fUB  = obj.computeFeasibleDesignVariable(lUB);                
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
        
        function fval = computeFeasibleDesignVariable(obj,lambda)
            obj.dualVariable.value = lambda;
            x                      = obj.updatePrimal();
            obj.updateDesignVariable(x);
            obj.constraint.computeFunctionAndGradient();
            fval = obj.constraint.value;
        end

        function x = updatePrimal(obj)
            lb = obj.lowerBound;
            ub = obj.upperBound;
            t  = obj.tau;
            Dg = obj.constraint.gradient;
            DJ = obj.cost.gradient;
            l  = obj.dualVariable.value;
            x  = obj.x;
            dx = -t*(DJ + l*Dg);
            xN = x + dx;
            x  = min(ub,max(xN,lb));
        end
        
        function updateDesignVariable(obj)
            obj.unconstrainedOptimizer.hasConverged = false;
            obj.unconstrainedOptimizer.compute();
            obj.unconstrainedOptimizer.updateConvergenceParams();
        end
        
    end
    
end