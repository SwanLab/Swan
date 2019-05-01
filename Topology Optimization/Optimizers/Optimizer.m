classdef Optimizer < handle
    
    properties (Access = public)
        convergenceVars
        targetParameters
    end
    
    properties (GetAccess = public, SetAccess = protected)
        historicalVariables
        niter
    end        
    
    properties (Access = protected)
        hasConverged
        hasFinished
        designVariable
        monitor
        cost
        constraint
        constraintCase
        historyPrinter

        maxiter 
        iStep
        nStep        
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        name
    end
    
    properties (Access = private)
        postProcess
        shallPrint
        printMode        
    end
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f = OptimizerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function solveProblem(obj,iStep,nStep)
            obj.iStep = iStep;
            obj.nStep = nStep;
            obj.createPostProcess();
            obj.cost.computeCostAndGradient();
            obj.constraint.computeCostAndGradient();
            obj.printOptimizerVariable();
            
            obj.hasFinished = false;
            
            while ~obj.hasFinished
                obj.niter = obj.niter+1;
                obj.update();
                obj.updateStatus();
                obj.refreshMonitoring()
                obj.printOptimizerVariable();
                obj.printHistory()
            end
            obj.hasConverged = 0;
        end
        
    end
    
    methods (Access = protected)
        
        function refreshMonitoring(obj)
           obj.monitor.refresh(obj.niter,obj.hasFinished,obj.iStep,obj.nStep);
        end
        
        function init(obj,cParams)
            obj.niter            = 0;
            obj.hasConverged     = false;            
            obj.constraintCase   = cParams.constraint_case;            
            obj.cost             = cParams.cost;
            obj.constraint       = cParams.constraint;            
            obj.designVariable   = cParams.designVar;
            obj.maxiter          = cParams.maxiter;
            obj.shallPrint       = cParams.shallPrint;
            obj.printMode        = cParams.printMode;
            
            obj.createHistoyPrinter();
            obj.createConvergenceVariables();
            obj.createMonitorDocker(cParams);
        end
        
      function updateStatus(obj)
            obj.hasFinished = obj.hasConverged || obj.niter >= obj.maxiter*(obj.iStep/obj.nStep);
        end
        
        function printOptimizerVariable(obj)
            d.x = obj.designVariable.value;
            d.cost = obj.cost;
            d.constraint = obj.constraint;
            obj.postProcess.print(obj.niter,d);
        end
        
        function printHistory(obj)
            obj.historyPrinter.print(obj.niter,obj.iStep);
        end
        
        function createPostProcess(obj)
            fileName = obj.designVariable.mesh.problemID;
            d = obj.createPostProcessDataBase(fileName);
            d.printMode  = obj.printMode;
            d.optimizer  = obj.name;
            d.cost       = obj.cost;
            d.constraint = obj.constraint;
            d.designVar  = obj.designVariable.type;
            if obj.shallPrint
                obj.postProcess = Postprocess('TopOptProblem',d);
            else
                obj.postProcess = Postprocess_Null('',d);
            end
        end        
        
    end
    
    methods (Access = private)
        
        function createHistoyPrinter(obj)
            cParams.fileName   = obj.designVariable.mesh.problemID;
            cParams.optimizer  = obj;
            cParams.cost       = obj.cost;
            cParams.constraint = obj.constraint;
            obj.historyPrinter = OptimizationMetricsPrinterFactory.create(obj.name,obj.shallPrint,cParams);
        end
        
        function createConvergenceVariables(obj)
            cVarD = ConvergenceVarsDispatcher.dispatchNames(obj.name);
            n = numel(cVarD);
            cVar = ConvergenceVariables(n);          
            obj.convergenceVars = cVar;            
        end
        
        function createMonitorDocker(obj,cParams)
            s = cParams.settingsMonitor;
            s.designVar       = obj.designVariable;
            s.optimizerName   = obj.name;
            s.cost            = obj.cost;
            s.constraint      = obj.constraint;
            s.convergenceVars = obj.convergenceVars;            
            obj.monitor = MonitoringDocker(s);
        end
        
        function d = createPostProcessDataBase(obj,fileName)
            d.mesh    = obj.designVariable.mesh;
            d.outName = fileName;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.getValue();
        end        
        
    end
    
end