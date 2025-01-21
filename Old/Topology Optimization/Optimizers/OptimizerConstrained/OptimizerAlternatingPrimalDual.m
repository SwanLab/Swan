classdef OptimizerAlternatingPrimalDual < Optimizer_PrimalDual
    
    properties (GetAccess = public, SetAccess = protected)
        type = 'AlternatingPrimalDual'
    end
    
    properties (Access = private)
        dualUpdater
    end
    
    methods (Access = public)
        
        function obj = OptimizerAlternatingPrimalDual(cParams)
            obj.init(cParams);
            obj.createLagrangian();
            obj.createOptimizerUnconstrained(cParams.uncOptimizerSettings);
            obj.createDualUpdater();
        end
        
        function update(obj)
            obj.lagrangian.updateBecauseOfPrimal();
            obj.updateDualVariable();
            obj.updatePrimalVariable();
            obj.updateConvergenceStatus();
        end
        
    end
    
    methods (Access = protected)
        
        function createLagrangianSettings(obj)
            cParams.type           = 'AugmentedLagrangian';
            cParams.constraintCase = obj.constraintCase;
            cParams.cost           = obj.cost;
            cParams.constraint     = obj.constraint;
            cParams.dualVariable   = obj.dualVariable;
            obj.lagrangianSettings = cParams;
        end
        
    end
    
    methods (Access = private)
        
        function createDualUpdater(obj)
            cParams.type                = 'AugmentedLagrangian';
            cParams.augmentedLagrangian = obj.lagrangian;
            cParams.constraint          = obj.constraint;
            cParams.cost                = obj.cost;
            cParams.dualVariable        = obj.dualVariable;
            obj.dualUpdater = DualUpdater.create(cParams);
        end
        
        function updateDualVariable(obj)
            obj.dualUpdater.updateDualVariable();
            obj.lagrangian.updateBecauseOfDual();
        end
        
        function updatePrimalVariable(obj)
            obj.unconstrainedOptimizer.update();
        end
        
    end
    
end