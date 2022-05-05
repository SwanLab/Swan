classdef Optimizer_PrimalDual < Optimizer
    
    properties (GetAccess = public, SetAccess = protected)
        unconstrainedOptimizer
    end
    
    properties (Access = protected)
        NSmerit
        NullSpaceSettings
        lagrangian
        lagrangianSettings
    end
    
    properties (Access = private)

    end
    
    methods (Access = public)
        
        function solveProblem(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
            obj.lagrangian.updateBecauseOfPrimal();
            obj.unconstrainedOptimizer.startLineSearch();
            obj.printOptimizerVariable();
            obj.hasFinished = false;

            while ~obj.hasFinished
                obj.increaseIter();
                obj.update();
                obj.updateStatus();
                obj.refreshMonitoring();
                obj.printOptimizerVariable();
                %obj.printHistory();
            end
            obj.printOptimizerVariable();
            obj.printHistory();

            obj.hasConverged = 0;
            obj.printHistoryFinalValues();
        end
        
    end
    
    methods (Access = protected)
        
        function updateOldValues(obj)
            obj.designVariable.updateOld();
            obj.dualVariable.updateOld();
            obj.cost.updateOld();
            obj.constraint.updateOld();
            obj.updateLagrangian();
            obj.lagrangian.updateOld();
        end
       
        function updateConvergenceStatus(obj)
            isOptimal   = obj.unconstrainedOptimizer.isOptimal();
            isFeasible  = obj.isFeasible();
            obj.hasConverged = isOptimal && isFeasible;
        end
        
        function itIs = isFeasible(obj)
            active_constr    = true(size(obj.dualVariable.value));
            constraintValues = abs(obj.constraint.value(active_constr));
            constrTol        = obj.targetParameters.constr_tol(active_constr);
            isNotFeasible = any(any(constraintValues > constrTol));
            itIs = ~isNotFeasible;
        end
        
        function createLagrangian(obj)
            obj.createLagrangianSettings();
            cParams = obj.lagrangianSettings;
            obj.lagrangian = ObjectiveFunction.create(cParams);
        end

        function createNullSpaceMerit(obj)
            obj.createNullSpaceSettings();
            cParams = obj.NullSpaceSettings();
            obj.NSmerit = ObjectiveFunction.create(cParams);
        end
        
       function createOptimizerUnconstrained(obj,cParams)
            cParams.lagrangian      = obj.lagrangian;
            cParams.convergenceVars = obj.convergenceVars;
            obj.unconstrainedOptimizer = Optimizer_Unconstrained.create(cParams);
       end 

       function createOptimizerUnconstrainedNS(obj,cParams)
           cParams.NSmerit         = obj.nsMerit;
           cParams.convergenceVars = obj.convergenceVars;
           obj.unconstrainedOptimizer = Optimizer_Unconstrained.create(cParams);
       end
       
        function updateLagrangian(obj)
            obj.lagrangian.computeFunction();
            obj.lagrangian.computeGradient();
        end
        
    end

    methods (Access = private)
       
        function increaseIter(obj)
            obj.nIter = obj.nIter+1;
        end
        
    end
    
    methods (Access = protected, Abstract)
        createLagrangianSettings(obj)
    end
    
end
