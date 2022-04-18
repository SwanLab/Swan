classdef Optimizer < handle
    
    properties (GetAccess = public, SetAccess = protected)
%         nIter
%         convergenceVars
%         historicalVariables
    end
    
    properties (Access = protected)
%         hasConverged
%         hasFinished
        designVariable
        dualVariable
        cost
        constraint
        outputFunction
        maxIter
        nIter
%         constraintCase
%         historyPrinter
%         targetParameters
%         incrementalScheme
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        type
    end
    
    properties (Access = private)
%         optimizerType
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f               = OptimizerFactory();
            obj             = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
%        function solveProblem(obj)
%             obj.cost.computeFunctionAndGradient();
%             obj.constraint.computeFunctionAndGradient();
%          %   obj.lagrangian.updateBecauseOfPrimal();
%          %   obj.unconstrainedOptimizer.startLineSearch();
%             obj.printOptimizerVariable();
%             obj.hasFinished = false;
% 
%             while ~obj.hasFinished
%                 obj.increaseIter();
%                 obj.update();
%                 obj.updateStatus();
%                 obj.refreshMonitoring();
%                 obj.printOptimizerVariable();
%                 %obj.printHistory();
%             end
%             obj.printOptimizerVariable();
%             obj.printHistory();
% 
%             obj.hasConverged = 0;
%             obj.printHistoryFinalValues();
%        end

       % PETARA PER TOT ARREU CAL ARREGLAR
        
    end
    
    methods (Access = protected)
        
        function initOptimizer(obj,cParams)
            obj.nIter             = 0;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.maxIter           = cParams.maxIter;
            obj.outputFunction    = cParams.outputFunction.monitoring;
        end
        
    end
    
end