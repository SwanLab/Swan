classdef OptimizerDualNestedInPrimal < Optimizer_PrimalDual
    
    properties (GetAccess = public, SetAccess = protected)
        name = 'DualNestedInPrimal'
    end
    
    properties (Access = private)
        desVarChangedValue
        costIncrease
        constraintProjector
    end
    
    methods (Access = public)
        
        function obj = OptimizerDualNestedInPrimal(cParams)
            obj.init(cParams);
            obj.createLagrangian();
            obj.createOptimizerUnconstrained(cParams.uncOptimizerSettings)
            obj.createConstraintProjector();
        end
        
        function solveProblem(obj)
            obj.cost.computeCostAndGradient();
            obj.constraint.computeCostAndGradient();
            obj.printOptimizerVariable();
            
            obj.designVariable.updateOld();
            obj.computeFeasibleDesignVariable();
                        
            obj.cost.computeCostAndGradient();

            obj.updateOldValues();
            
            obj.unconstrainedOptimizer.updateConvergenceParams();
            obj.refreshMonitoring();
            obj.printOptimizerVariable();
            obj.printHistory();
            
            obj.hasFinished = false;
            
            while ~obj.hasFinished
                obj.niter = obj.niter+1;
                
                obj.unconstrainedOptimizer.initLineSearch();
                kappa = 0.1*obj.unconstrainedOptimizer.line_search.kappa;
                obj.unconstrainedOptimizer.line_search.kappa = kappa;
                
                while ~obj.hasUnconstraintedOptimizerConverged()
                    obj.restartValues();
                    obj.computeValue();
                    obj.unconstrainedOptimizer.line_search.computeKappa();
                end
                
                obj.unconstrainedOptimizer.updateConvergenceParams();
                
                obj.updateConvergenceStatus();
                obj.updateStatus();
                
                obj.updateOldValues();
                
                obj.refreshMonitoring();
                obj.printOptimizerVariable();
                obj.printHistory();
            end
            obj.hasConverged = 0;
        end
        
        function restartValues(obj)
            obj.designVariable.restart();
            obj.dualVariable.restart();
            obj.cost.restart();
            obj.constraint.restart();
            obj.updateLagrangian();
        end
        
        function updateOldValues(obj)
            obj.designVariable.updateOld();
            obj.dualVariable.updateOld();
            obj.cost.updateOld();
            obj.constraint.updateOld();
            obj.updateLagrangian();
            obj.lagrangian.updateOld();
        end
        
        function computeFeasibleDesignVariable(obj)
            if obj.isNotFeasible()
                obj.constraint.updateOld();
                obj.updateLagrangian();
                obj.lagrangian.updateOld();
                
                %kappa
                obj.unconstrainedOptimizer.initLineSearch();
                kappa = 0.01*obj.unconstrainedOptimizer.line_search.kappa;
                obj.unconstrainedOptimizer.line_search.kappa = kappa;
                
                obj.constraintProjector.project();
            end
        end
        
        function itIsNot = isNotFeasible(obj)
            itIsNot = ~obj.isFeasible();
        end
        
        
        function update(obj)
            while ~obj.hasUnconstraintedOptimizerConverged()
                obj.computeValue();
                obj.unconstrainedOptimizer.line_search.computeKappa();
            end
        end
        
    end
    
    methods (Access = protected)
        
        function createLagrangianSettings(obj)
            cParams.type         = 'Lagrangian';
            cParams.cost         = obj.cost;
            cParams.constraint   = obj.constraint;
            cParams.dualVariable = obj.dualVariable;
            obj.lagrangianSettings = cParams;
        end
        
    end
    
    methods (Access = private)
        
        function createConstraintProjector(obj)
            cParams.cost           = obj.cost;
            cParams.constraint     = obj.constraint;
            cParams.designVariable = obj.designVariable;
            cParams.dualVariable   = obj.dualVariable;
            cParams.lagrangian     = obj.lagrangian;
            cParams.targetParameters = obj.targetParameters;
            cParams.unconstrainedOptimizer = obj.unconstrainedOptimizer;
            obj.constraintProjector = ConstraintProjector(cParams);
        end
        
        function computeValue(obj)
            obj.constraintProjector.project();
            obj.cost.computeCostAndGradient();
            obj.updateLagrangian();
        end
        
        function itHas = hasUnconstraintedOptimizerConverged(obj)
            itHas = obj.isStepAcceptable() || obj.isLineSeachTooSmall();
        end
        
        function itIs = isStepAcceptable(obj)
            incr = obj.lagrangian.computeIncrement();
            costHasDecreased = incr < 0;
            itIs = costHasDecreased;
        end
        
        function itIs = isLineSeachTooSmall(obj)
            kappa     = obj.unconstrainedOptimizer.line_search.kappa;
            kappa_min = obj.unconstrainedOptimizer.line_search.kappa_min;
            itIs = kappa <= kappa_min;
        end
        
        function updateLagrangian(obj)
            obj.lagrangian.computeFunction();
            obj.lagrangian.computeGradient();
        end
        
    end
end