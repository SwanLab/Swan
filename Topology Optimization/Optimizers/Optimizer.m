classdef Optimizer < handle

    properties (GetAccess = public, SetAccess = protected)
        niter
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

        maxiter     
        incrementalScheme                
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        name
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
            obj.cost.computeCostAndGradient();
            obj.constraint.computeCostAndGradient();
            obj.printOptimizerVariable();
            
            obj.hasFinished = false;
            
            while ~obj.hasFinished
                obj.niter = obj.niter+1;
                obj.update();
                obj.updateStatus();
                obj.refreshMonitoring();
                obj.printOptimizerVariable();
                obj.printHistory();
            end
            obj.hasConverged = 0;
        end
        
    end
    
    methods (Access = protected)
        
        function refreshMonitoring(obj)
           iStep = obj.incrementalScheme.iStep;
           nStep = obj.incrementalScheme.nSteps;
           obj.monitor.refresh(obj.niter,obj.hasFinished,iStep,nStep);
        end
        
        function init(obj,cParams)
            obj.niter            = 0;
            obj.hasConverged     = false;            
            obj.constraintCase   = cParams.constraint_case;            
            obj.cost             = cParams.cost;
            obj.constraint       = cParams.constraint;            
            obj.designVariable   = cParams.designVar;
            obj.dualVariable     = cParams.dualVariable;
            obj.maxiter          = cParams.maxiter;
            obj.incrementalScheme = cParams.incrementalScheme;
            obj.targetParameters = cParams.target_parameters;
            
            obj.createHistoryPrinter(cParams.historyPrinterSettings);
            obj.createConvergenceVariables();
            obj.createMonitorDocker(cParams.settingsMonitor);
            obj.createPostProcess(cParams.postProcessSettings);            
        end
        
      function updateStatus(obj)
            iStep = obj.incrementalScheme.iStep;
            nStep = obj.incrementalScheme.nSteps;
            obj.hasFinished = obj.hasConverged || obj.niter >= obj.maxiter*(iStep/nStep);
        end
        
        function printOptimizerVariable(obj)
            d.x = obj.designVariable.value;
            d.cost = obj.cost;
            d.constraint = obj.constraint;
            obj.postProcess.print(obj.niter,d);
        end
        
        function printHistory(obj)
           iStep = obj.incrementalScheme.iStep;
           obj.historyPrinter.print(obj.niter,iStep);
        end
        
    end
    
    methods (Access = private)
        
        function createHistoryPrinter(obj,cParams)
            cParams.optimizer  = obj;
            cParams.cost       = obj.cost;
            cParams.constraint = obj.constraint;
            obj.historyPrinter = OptimizationMetricsPrinterFactory.create(cParams);
        end
        
        function createConvergenceVariables(obj)
            cVarD = ConvergenceVarsDispatcher.dispatchNames(obj.name);
            n = numel(cVarD);
            cVar = ConvergenceVariables(n);          
            obj.convergenceVars = cVar;            
        end
        
        function createMonitorDocker(obj,s)
            s.designVar       = obj.designVariable;
            s.cost            = obj.cost;
            s.constraint      = obj.constraint;
            s.convergenceVars = obj.convergenceVars;            
            obj.monitor = MonitoringDocker(s);
        end
        
        function createPostProcess(obj,cParams)
            d = obj.createPostProcessDataBase();
            d.printMode = cParams.printMode;
            if cParams.shallPrint
                obj.postProcess = Postprocess('TopOptProblem',d);
            else
                obj.postProcess = Postprocess_Null('',d);
            end
        end        
        
        function d = createPostProcessDataBase(obj)
            d.mesh    = obj.designVariable.mesh;
            d.outName = obj.designVariable.mesh.problemID;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.getValue();
            d.optimizer  = obj.name;
            d.cost       = obj.cost;
            d.constraint = obj.constraint;
            d.designVar  = obj.designVariable.type;            
        end        
        
    end
    
end