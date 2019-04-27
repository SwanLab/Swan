classdef Optimizer_Constrained < Optimizer
    
    properties (Access = public)
        fhtri
        niter = 0
        maxiter
    end
    
    properties (Access = protected)
        hasFinished
        designVar
        monitor
        cost
        constraint
        constraintCase
        historyPrinter
        iStep
        nStep
    end
    
    properties (Access = private)
        postProcess
        printing
        printMode
    end
    
    methods (Access = public)
        
        function solveProblem(obj,iStep,nStep)
            obj.iStep = iStep;
            obj.nStep = nStep;
            obj.createPostProcess();
            x0 = obj.designVar.value;
            obj.cost.computeCostAndGradient(x0);
            obj.constraint.computeCostAndGradient(x0);
            obj.printOptimizerVariable();
            
            obj.hasFinished = false;
            
             while ~obj.hasFinished
                obj.niter = obj.niter+1;
                obj.update();
                obj.updateStatus();
                obj.monitor.refresh(obj.niter,obj.hasFinished,iStep,nStep);
                obj.printOptimizerVariable();
                obj.printHistory()
            end
            obj.printOptimizerVariable();            
            obj.hasConverged = 0;
        end
        
    end
    
    methods (Access = protected)
        
        function updateStatus(obj)
            obj.hasFinished = obj.hasConverged || obj.niter >= obj.maxiter*(obj.iStep/obj.nStep);
        end
        
        function printOptimizerVariable(obj)
            d.x = obj.designVar.value;
            d.cost = obj.cost;
            d.constraint = obj.constraint;
            obj.postProcess.print(obj.niter,d);
        end
        
        function printHistory(obj)
            obj.historyPrinter.print(obj.niter,obj.iStep);
        end
        
        function createPostProcess(obj)
            fileName = obj.designVar.mesh.problemID;
            d = obj.createPostProcessDataBase(fileName);
            d.printMode = obj.printMode;
            d.optimizer = obj.name;
            d.cost = obj.cost;
            d.constraint = obj.constraint;
            if obj.printing
                obj.postProcess = Postprocess('TopOptProblem',d);
            else
                obj.postProcess = Postprocess_Null('',d);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            set = cParams.settings;
            
            obj.constraintCase   = set.constraint_case;
            obj.hasConverged     = false;
            
            obj.cost = cParams.cost;
            obj.constraint = cParams.constraint;
            
            obj.designVar = cParams.designVar;
            obj.maxiter   = set.maxiter;
            obj.printing  = set.printing;
            obj.printMode = set.printMode;
            obj.createHistoyPrinter();
            
            cParams.monitorSettings.showOptParams               = cParams.monitoring;
            cParams.monitorSettings.refreshInterval             = cParams.monitoring_interval;
            cParams.monitorSettings.shallDisplayDesignVar       = cParams.plotting;
            cParams.monitorSettings.shallShowBoundaryConditions = cParams.showBC;
            
            cParams.monitorSettings.problemID                   = set.case_file;
            cParams.monitorSettings.costFuncNames               = set.cost;
            cParams.monitorSettings.costWeights                 = set.weights;
            cParams.monitorSettings.constraintFuncs             = set.constraint;
            cParams.monitorSettings.dim                         = set.pdim;
            
            obj.createMonitorDocker(cParams.monitorSettings);
        end
        
    end
    
    methods (Access = private)
        
        function createHistoyPrinter(obj)
            settingsMetricsPrinter.fileName = obj.designVar.mesh.problemID;
            settingsMetricsPrinter.optimizer = obj;
            settingsMetricsPrinter.cost = obj.cost;
            settingsMetricsPrinter.constraint = obj.constraint;
            obj.historyPrinter = OptimizationMetricsPrinterFactory.create(obj.name,obj.printing,settingsMetricsPrinter);
        end
        
        function createMonitorDocker(obj,cParams)
            mS.settingsParamsMonitor.showOptParams    = cParams.showOptParams;
            mS.settingsParamsMonitor.refreshInterval  = cParams.refreshInterval;
            mS.settingsParamsMonitor.problemID        = cParams.problemID;
            mS.settingsParamsMonitor.costFuncNames    = cParams.costFuncNames;
            mS.settingsParamsMonitor.costWeights      = cParams.costWeights;
            mS.settingsParamsMonitor.constraintFuncs  = cParams.constraintFuncs;
            mS.settingsParamsMonitor.optimizerName    = obj.name;
            mS.settingsParamsMonitor.cost             = obj.cost;
            mS.settingsParamsMonitor.constraint       = obj.constraint;
            mS.settingsParamsMonitor.convergenceVars  = obj.convergenceVars;

            mS.settingsDesignVarMonitor.shallDisplay  = cParams.shallDisplayDesignVar;
            mS.settingsDesignVarMonitor.showBC        = cParams.shallShowBoundaryConditions;
            mS.settingsDesignVarMonitor.designVar     = obj.designVar;
            mS.settingsDesignVarMonitor.optimizerName = obj.name;
            mS.settingsDesignVarMonitor.dim           = cParams.dim;
            
            obj.monitor = MonitoringDocker(mS);
        end
        
        function d = createPostProcessDataBase(obj,fileName)
            d.mesh    = obj.designVar.mesh;
            d.outName = fileName;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.getValue();
        end
        
    end
    
end
