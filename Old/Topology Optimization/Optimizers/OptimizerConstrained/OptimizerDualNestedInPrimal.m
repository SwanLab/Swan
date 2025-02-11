classdef OptimizerDualNestedInPrimal < Optimizer_PrimalDual

    properties (GetAccess = public, SetAccess = protected)
        type = 'DualNestedInPrimal'
    end

    properties (Access = private)
        desVarChangedValue
        costIncrease
        constraintProjector
    end

    methods (Access = public)

        function obj = OptimizerDualNestedInPrimal(cParams)
            obj.initOptimizer(cParams);
            obj.createLagrangian();
            obj.createOptimizerUnconstrained(cParams.uncOptimizerSettings)
            obj.createConstraintProjector();
        end

        function solveProblem(obj)
            obj.hasFinished = false;
            
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            
            
            obj.updateOldValues();
            obj.unconstrainedOptimizer.startLineSearch();
            obj.unconstrainedOptimizer.updateConvergenceParams();
            obj.refreshMonitoring();
            obj.printOptimizerVariable();
       %     obj.printHistory();
            obj.nIter = obj.nIter+1;
            
            
            obj.designVariable.updateOld();
            obj.computeFeasibleDesignVariable();

            obj.cost.computeFunctionAndGradient();

            obj.updateOldValues();

            obj.unconstrainedOptimizer.updateConvergenceParams();   
            obj.unconstrainedOptimizer.updateLineSearch();
            obj.refreshMonitoring();
            obj.printHistory();
            obj.saveDesignVariable();
            obj.printOptimizerVariable();
       %     obj.printHistory();
            obj.nIter = obj.nIter+1;

            %obj.hasFinished = false;
            obj.updateStatus();

            
            while ~obj.hasFinished
             

                obj.unconstrainedOptimizer.tryLineSearch();
                while ~obj.hasUnconstraintedOptimizerConverged()
                    obj.restartValues();
                    obj.computeValue();
                    if ~obj.hasUnconstraintedOptimizerConverged()
                        obj.unconstrainedOptimizer.updateLineSearch();
                    end
                end

                obj.cost.computeFunctionAndGradient();

                obj.unconstrainedOptimizer.updateConvergenceParams();

                obj.updateConvergenceStatus();
                obj.updateStatus();

                obj.updateOldValues();

                obj.refreshMonitoring();
                obj.printOptimizerVariable();
                obj.printHistory();
                obj.saveDesignVariable();
                obj.nIter = obj.nIter+1;
            end
            %obj.printOptimizerVariable();
            %obj.printHistory();
            obj.hasConverged = 0;
            obj.printHistoryFinalValues();
        end
        
        function saveDesignVariable(obj)
            x = obj.designVariable.value;
            mesh = obj.designVariable.mesh.innerMeshOLD;
            path = 'Output/CantileverTetraPerimeterTotal/DesignVariable';
            %save([path,num2str(obj.nIter)],'x','mesh');
        end

        function restartValues(obj)
            obj.designVariable.restart();
            obj.dualVariable.restart();
            obj.cost.restart();
            obj.constraint.restart();
            obj.updateLagrangian();
        end



        function computeFeasibleDesignVariable(obj)
            if obj.isNotFeasible()
                obj.constraint.updateOld();
                obj.updateLagrangian();
                obj.lagrangian.updateOld();

                obj.unconstrainedOptimizer.tryLineSearch();
                %kappa = obj.unconstrainedOptimizer.lineSearch.kappa;
                %obj.unconstrainedOptimizer.lineSearch.kappa = kappa;

                obj.constraintProjector.project();
            end
        end

        function itIsNot = isNotFeasible(obj)
            itIsNot = ~obj.isFeasible();
        end


        function update(obj)
            while ~obj.hasUnconstraintedOptimizerConverged()
                obj.computeValue();
                obj.unconstrainedOptimizer.updateLineSearch();
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
        
        function updateLagrangian(obj)
            obj.lagrangian.computeFunction();
            obj.lagrangian.computeGradient();
        end
        
        function updateOldValues(obj)
            obj.designVariable.updateOld();
            obj.dualVariable.updateOld();
            obj.cost.updateOld();
            obj.constraint.updateOld();
            obj.updateLagrangian();
            obj.lagrangian.updateOld();
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
            obj.cost.computeFunction();
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
            itIs = obj.unconstrainedOptimizer.isLineSearchTooSmall();
        end
             


    end
end
