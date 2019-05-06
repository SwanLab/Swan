classdef Optimizer_AugLag < Optimizer_PrimalDual
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'AugmentedLagrangian'
    end
    
    properties (GetAccess = public, SetAccess = private)
        unconstrainedOptimizer
    end
    
    properties (Access = private)
        augLagrangian      
        dualUpdater
    end
    
    methods (Access = public)
        
        function obj = Optimizer_AugLag(cParams)
            obj.init(cParams);            
            obj.createAugmentedLagrangian();
            cParams.uncOptimizerSettings.lagrangian = obj.augLagrangian;           
            cParams.uncOptimizerSettings.convergenceVars = obj.convergenceVars;
            obj.unconstrainedOptimizer = Optimizer_Unconstrained.create(cParams.uncOptimizerSettings);                           
            obj.createDualUpdater()
        end
        
        function createDualUpdater(obj)
            cParams.type                = 'AugmentedLagrangian';
            cParams.augmentedLagrangian = obj.augLagrangian;
            cParams.constraint          = obj.constraint;
            cParams.dualVariable        = obj.dualVariable;
            obj.dualUpdater = DualUpdater.create(cParams);
        end
        
        function update(obj)
            obj.updateDualVariable();
            obj.augLagrangian.updateBecauseOfDual();
            obj.updatePrimalVariable();
            obj.updateConvergenceStatus();
        end
        
    end
    
    methods (Access = private)
        
        function updatePrimalVariable(obj)
            obj.designVariable.updateOld();                        
            obj.unconstrainedOptimizer.init();            
            while ~obj.unconstrainedOptimizer.hasConverged
                obj.designVariable.restart();
                obj.unconstrainedOptimizer.update();
            end
            obj.revertIfDesignNotImproved();
        end
        
        function updateConvergenceStatus(obj)
            active_constr = true(size(obj.dualVariable.value));
            isNotOptimal  = obj.unconstrainedOptimizer.opt_cond >=  obj.unconstrainedOptimizer.optimality_tol;
            isNotFeasible = any(any(abs(obj.constraint.value(active_constr)) > obj.unconstrainedOptimizer.constr_tol(active_constr)));
            hasNotConverged = isNotOptimal || isNotFeasible;
            obj.hasConverged = ~hasNotConverged;
        end
        
        function updateDualVariable(obj)            
            obj.dualUpdater.updateDualVariable();
        end
        
        function revertIfDesignNotImproved(obj)
            if ~obj.unconstrainedOptimizer.designImproved
                obj.designVariable.restart();
            end
        end
        
        function createAugmentedLagrangian(obj)
            cParams.constraintCase = obj.constraintCase;
            cParams.cost           = obj.cost;
            cParams.constraint     = obj.constraint;
            cParams.dualVariable   = obj.dualVariable;
            obj.augLagrangian = AugmentedLagrangian(cParams);
        end

    end
    
end