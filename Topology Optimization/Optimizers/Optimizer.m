classdef Optimizer < handle
    
    properties (GetAccess = public, SetAccess = protected)
        nIter
        convergenceVars
        historicalVariables
    end
    
    properties (Access = protected)
        hasConverged
        hasFinished
        designVariable
        dualVariable
        cost
        constraint
        constraintCase
        historyPrinter
        targetParameters
        
        maxIter
        incrementalScheme
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        type
    end
    
    properties (Access = private)
        postProcess
        monitor
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = OptimizerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
       function solveProblem(obj)
            obj.cost.computeFunctionAndGradient();
            obj.constraint.computeFunctionAndGradient();
         %   obj.lagrangian.updateBecauseOfPrimal();            
         %   obj.unconstrainedOptimizer.startLineSearch();            
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
        
        function saveMonitoring(obj)
            obj.monitor.saveMonitorFigure();
        end         
        
    end
    
    methods (Access = protected)
        
        function refreshMonitoring(obj)
            iStep = obj.incrementalScheme.iStep;
            nStep = obj.incrementalScheme.nSteps;
            obj.monitor.refresh(obj.nIter,obj.hasFinished,iStep,nStep);
        end
        
        function init(obj,cParams)
            obj.nIter             = 0;
            obj.hasConverged      = false;
            obj.constraintCase    = cParams.constraintCase;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.maxIter           = cParams.maxIter;
            obj.incrementalScheme = cParams.incrementalScheme;
            obj.targetParameters  = cParams.targetParameters;
            
            obj.createHistoryPrinter(cParams.historyPrinterSettings);
            obj.createConvergenceVariables(cParams);
            obj.createMonitorDocker(cParams.monitoringDockerSettings);
            obj.createPostProcess(cParams.postProcessSettings);
        end
        
        function updateStatus(obj)
            obj.hasFinished = obj.hasConverged || obj.hasExceededStepIterations();
        end
        
        function printOptimizerVariable(obj)
            d.fields  = obj.designVariable.getVariablesToPlot();
            d.cost = obj.cost;
            d.constraint = obj.constraint;
            obj.postProcess.print(obj.nIter,d);
        end
        
        function printHistory(obj)
            iStep = obj.incrementalScheme.iStep;
            obj.historyPrinter.print(obj.nIter,iStep);
        end
        
        function printHistoryFinalValues(obj)
            obj.historyPrinter.printFinal();                    
        end
               
    end
    
    methods (Access = private)
        
        function increaseIter(obj)
            obj.nIter = obj.nIter+1;
        end 
        
        function itHas = hasExceededStepIterations(obj)
            iStep = obj.incrementalScheme.iStep;
            nStep = obj.incrementalScheme.nSteps;
            itHas = obj.nIter >= obj.maxIter*(iStep/nStep);
        end
        
        function createHistoryPrinter(obj,cParams)
            cParams.optimizer  = obj;
            cParams.cost       = obj.cost;
            cParams.constraint = obj.constraint;
            obj.historyPrinter = OptimizationMetricsPrinterFactory.create(cParams);
        end
        
        function createConvergenceVariables(obj,cParams)
            s = cParams.optimizerNames;
            cVarD = ConvergenceVarsDispatcher.dispatchNames(s);
            n = numel(cVarD);
            cVar = ConvergenceVariables(n);
            obj.convergenceVars = cVar;
        end
        
        function createMonitorDocker(obj,s)
            s.designVariable  = obj.designVariable;
            s.dualVariable    = obj.dualVariable;
            s.cost            = obj.cost;
            s.constraint      = obj.constraint;
            s.convergenceVars = obj.convergenceVars;
            obj.monitor = MonitoringDocker(s);
        end
        
        function createPostProcess(obj,cParams)
            d = obj.createPostProcessDataBase(cParams);
            d.printMode = cParams.printMode;
            d.nDesignVariables = obj.designVariable.nVariables;
            if cParams.shallPrint
                obj.postProcess = Postprocess('TopOptProblem',d);
            else
                obj.postProcess = Postprocess_Null('',d);
            end
        end
        
        function d = createPostProcessDataBase(obj,cParams)
            d.mesh    = obj.designVariable.mesh;
            d.outName = cParams.femFileName;
            d.pdim    = cParams.pdim;
            d.ptype   = cParams.ptype;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.getValue();
            d.optimizer  = obj.type;
            d.cost       = obj.cost;
            d.constraint = obj.constraint;
            d.designVar  = obj.designVariable.type;
        end
        

        
    end
    
end