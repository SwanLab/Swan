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
            tolCons = obj.targetParameters.constr_tol;
            obj.problem.options = optimset(obj.problem.options,'TolX',1e-2*tolCons);            
            fzero(obj.problem);
        end
        
    end
    
    methods (Access = private)
        
        function defineProblem(obj)
            obj.problem.solver = 'fzero';
            obj.problem.options = optimset(@fzero);
            obj.problem.objective = @(lambda) obj.computeFeasibleDesignVariable(lambda);
            obj.problem.x0 = [0 100];
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
        
        function restartValues(obj)

        end
        
    end
 
end