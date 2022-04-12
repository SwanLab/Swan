classdef Optimizer < handle
    
    properties (GetAccess = public, SetAccess = protected)
%         nIter
%         convergenceVars
%         historicalVariables
    end
    
    properties (Access = protected)
%         hasConverged
%         hasFinished
%         designVariable
%         dualVariable
%         cost
%         constraint
%         constraintCase
%         historyPrinter
%         targetParameters
%         outputFunction
%         maxIter
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
        
%         function init(obj,cParams)
%             obj.nIter             = 0;
%             obj.hasConverged      = false;
%             obj.constraintCase    = cParams.constraintCase;
%             obj.cost              = cParams.cost;
%             obj.constraint        = cParams.constraint;
%             obj.designVariable    = cParams.designVar;
%             obj.dualVariable      = cParams.dualVariable;
%             obj.maxIter           = cParams.maxIter;
%             obj.incrementalScheme = cParams.incrementalScheme;
%             obj.targetParameters  = cParams.targetParameters;
%             obj.optimizerType     = cParams.outputFunction.type;
%             obj.outputFunction.monitoring.create(cParams);
%         end

%         function updateIterInfo(obj)
%             switch obj.optimizerType
%                 case 'Topology'
%                     obj.increaseIter();
%                     obj.updateStatus();
%                 otherwise
%             end
%         end
        
    end
    
    methods (Access = private)

%         function increaseIter(obj)
%             obj.nIter = obj.nIter+1;
%         end
% 
%         function updateStatus(obj)
%             obj.hasFinished = obj.hasConverged || obj.hasExceededStepIterations();
%         end
%         
%         function itHas = hasExceededStepIterations(obj)
%             iStep = obj.incrementalScheme.iStep;
%             nStep = obj.incrementalScheme.nSteps;
%             itHas = obj.nIter >= obj.maxIter*(iStep/nStep);
%         end
        
%         function d = createPostProcessDataBase(obj,cParams)
%             d.mesh    = obj.designVariable.mesh;
%             d.outName = cParams.femFileName;
%             d.pdim    = cParams.pdim;
%             d.ptype   = cParams.ptype;
%             ps = PostProcessDataBaseCreator(d);
%             d  = ps.getValue();
%             d.optimizer  = obj.type;
%             d.cost       = obj.cost;
%             d.constraint = obj.constraint;
%             d.designVar  = obj.designVariable.type;
%         end
        
    end
    
end