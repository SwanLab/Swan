classdef OptimizerAlternatingPrimalDual < Optimizer_PrimalDual
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'AlternatingPrimalDual'
    end
    
    properties (GetAccess = public, SetAccess = private)
        unconstrainedOptimizer
    end
    
    properties (Access = private)
        augLagrangian      
        dualUpdater
    end
    
    methods (Access = public)
        
        function obj = OptimizerAlternatingPrimalDual(cParams)
            obj.init(cParams);            
            obj.createAugmentedLagrangian();
            obj.createOptimizerUnconstrained(cParams.uncOptimizerSettings);                        
            obj.createDualUpdater();
        end
        
        function update(obj)
            obj.updateDualVariable();
            obj.updatePrimalVariable();
            obj.updateConvergenceStatus();
        end
        
    end
    
    methods (Access = private)
        
        function createAugmentedLagrangian(obj)
            cParams.constraintCase = obj.constraintCase;
            cParams.cost           = obj.cost;
            cParams.constraint     = obj.constraint;
            cParams.dualVariable   = obj.dualVariable;
            obj.augLagrangian = AugmentedLagrangian(cParams);
        end        
        
        function createOptimizerUnconstrained(obj,cParams)
            cParams.lagrangian      = obj.augLagrangian;           
            cParams.convergenceVars = obj.convergenceVars;
            obj.unconstrainedOptimizer = Optimizer_Unconstrained.create(cParams);              
        end       
        
        function createDualUpdater(obj)
            cParams.type                = 'AugmentedLagrangian';
            cParams.augmentedLagrangian = obj.augLagrangian;
            cParams.constraint          = obj.constraint;
            cParams.dualVariable        = obj.dualVariable;
            obj.dualUpdater = DualUpdater.create(cParams);
        end   
        
        function updateDualVariable(obj)            
            obj.dualUpdater.updateDualVariable();
            obj.augLagrangian.updateBecauseOfDual();            
        end        
        
        function updatePrimalVariable(obj)
            obj.unconstrainedOptimizer.update();
        end
        
        function updateConvergenceStatus(obj)
            active_constr = true(size(obj.dualVariable.value));
            isNotOptimal  = obj.unconstrainedOptimizer.opt_cond >=  obj.unconstrainedOptimizer.optimality_tol;
            isNotFeasible = any(any(abs(obj.constraint.value(active_constr)) > obj.unconstrainedOptimizer.constr_tol(active_constr)));
            hasNotConverged = isNotOptimal || isNotFeasible;
            obj.hasConverged = ~hasNotConverged;
        end
        
    end
    
end