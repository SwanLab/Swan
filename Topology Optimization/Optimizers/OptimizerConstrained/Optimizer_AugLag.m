classdef Optimizer_AugLag < Optimizer_Constrained
    
    properties (GetAccess = public, SetAccess = private)
        unconstrainedOptimizer
        augLagrangian
        penalty
    end
    
    properties (Access = private)
        x
        lambda
    end
    
    methods (Access = public)
        
        function obj = Optimizer_AugLag(cParams)
            obj.init(cParams);
            
            obj.unconstrainedOptimizer = cParams.unconstrainedOptimizer;
            obj.convergenceVars = obj.unconstrainedOptimizer.convergenceVars;
            
            obj.createAugmentedLagrangian();
            obj.createLambdaAndPenalty();
        end
        
        function x = update(obj,x0)
            obj.unconstrainedOptimizer.target_parameters = obj.target_parameters;
            obj.updateDualVariable();
            obj.augLagrangian.updateBecauseOfDual(obj.lambda,obj.penalty);
            obj.updatePrimalVariable(x0);
            obj.updateConvergenceStatus();
            x = obj.x;
        end
        
    end
    
    methods (Access = private)
        
        function updatePrimalVariable(obj,x0)
            obj.unconstrainedOptimizer.init(x0,obj.augLagrangian);
            while ~obj.unconstrainedOptimizer.hasConverged
                x = obj.unconstrainedOptimizer.update(x0);
            end
            x = obj.revertIfDesignNotImproved(x,x0);
            
            obj.x = x;
        end
        
        function updateConvergenceStatus(obj)
            active_constr = obj.penalty > 0;
            isNotOptimal  = obj.unconstrainedOptimizer.opt_cond >=  obj.unconstrainedOptimizer.optimality_tol;
            isNotFeasible = any(any(abs(obj.constraint.value(active_constr)) > obj.unconstrainedOptimizer.constr_tol(active_constr)));
            hasNotConverged = isNotOptimal || isNotFeasible;
            obj.hasConverged = ~hasNotConverged;
        end
        
        function updateDualVariable(obj)
            l   = obj.lambda;
            rho = obj.penalty;
            c   = obj.constraint.value';
            l = l + rho.*c;
            obj.lambda = l;
        end
        
        function x = revertIfDesignNotImproved(obj,x,x0)
            if ~obj.unconstrainedOptimizer.designImproved
                x = x0;
            end
        end
        
        function createAugmentedLagrangian(obj)
            augLagS.constraintCase = obj.constraintCase;
            obj.augLagrangian = AugmentedLagrangian(augLagS);
            obj.augLagrangian.link(obj.cost,obj.constraint);
        end
        
        function createLambdaAndPenalty(obj)
            nConstraints = obj.constraint.nSF;
            obj.lambda  = zeros(1,nConstraints);
            obj.penalty = ones(1,nConstraints);
        end
        
    end
    
end