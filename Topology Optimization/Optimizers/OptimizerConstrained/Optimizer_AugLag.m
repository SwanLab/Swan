classdef Optimizer_AugLag < Optimizer_Constrained
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'AugmentedLagrangian'
    end
    
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
            obj.createAugmentedLagrangian();
            obj.createLambdaAndPenalty();
            obj.unconstrainedOptimizer = Optimizer_Unconstrained.create(cParams);                           
            obj.unconstrainedOptimizer.init(obj.augLagrangian);            
        end
        
        function update(obj)
            obj.updateDualVariable();
            obj.augLagrangian.updateBecauseOfDual(obj.lambda,obj.penalty);
            obj.updatePrimalVariable(obj.designVariable.value);
            obj.updateConvergenceStatus();
        end
        
    end
    
    methods (Access = private)
        
        function updatePrimalVariable(obj,x0)
            obj.designVariable.valueOld = x0;
            obj.unconstrainedOptimizer.init2();            
            while ~obj.unconstrainedOptimizer.hasConverged
                obj.designVariable.value = obj.designVariable.valueOld;
                obj.unconstrainedOptimizer.update();
            end
            obj.revertIfDesignNotImproved();
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
        
        function revertIfDesignNotImproved(obj)
            if ~obj.unconstrainedOptimizer.designImproved
                obj.designVariable.value = obj.designVariable.valueOld;
            end
        end
        
        function createAugmentedLagrangian(obj)
            cParamsAL.constraintCase = obj.constraintCase;
            cParamsAL.cost           = obj.cost;
            cParamsAL.constraint     = obj.constraint;
            obj.augLagrangian = AugmentedLagrangian(cParamsAL);
        end
        
        function createLambdaAndPenalty(obj)
            nConstraints = obj.constraint.nSF;
            obj.lambda  = zeros(1,nConstraints);
            obj.penalty = ones(1,nConstraints);
        end
        
    end
    
end