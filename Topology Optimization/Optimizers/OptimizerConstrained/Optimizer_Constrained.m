classdef Optimizer_Constrained < Optimizer
    
    properties (Access = public)
        fhtri
        niter = 0
        optimizer
        maxiter
    end
    
    properties (Access = protected)
        designVar
        monitor
        cost
        constraint
        constraintCase
        historyPrinter
    end
    
    properties (Access = private)
        postProcess
        printing
        printMode
    end
    
    methods (Access = public)
        
        function designVar = solveProblem(obj,designVar,istep,nstep)
            obj.createPostProcess();
            x0 = designVar.value;
            obj.cost.computeCostAndGradient(x0);
            obj.constraint.computeCostAndGradient(x0);
            obj.printOptimizerVariable();
            
            %obj.monitor.refresh(x0,obj.niter,obj.cost,obj.constraint,obj.hasFinished(istep,nstep),istep,nstep);
            
            while ~obj.hasFinished(istep,nstep)
                obj.niter = obj.niter+1;
                x = obj.update(x0);
                designVar.update(x);
                obj.monitor.refresh(obj.niter,obj.cost,obj.constraint,obj.hasFinished(istep,nstep),istep,nstep);
                obj.printOptimizerVariable();
                obj.printHistory(istep)
                x0 = x;
            end
            obj.printOptimizerVariable();
            
            obj.hasConverged = 0;
        end
        
    end
    
    methods (Access = protected)
        
        function printOptimizerVariable(obj)
            d.x = obj.designVar.value;
            d.cost = obj.cost;
            d.constraint = obj.constraint;
            obj.postProcess.print(obj.niter,d);
        end
        
        function printHistory(obj,iStep)
            obj.historyPrinter.print(obj.niter,iStep);
        end
        
        function createPostProcess(obj)
            fileName = obj.designVar.meshGiD.problemID;
            d = obj.createPostProcessDataBase(fileName);
            d.printMode = obj.printMode;
            d.optimizer = obj.optimizer;
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
            set = cParams.settings.settings;
            designVar = cParams.designVariable;
            
            obj.constraintCase   = set.constraint_case;
            obj.hasConverged     = false;
            
            obj.cost = cParams.settings.cost;
            obj.constraint = cParams.settings.constraint;
            
            obj.designVar = designVar;
            obj.optimizer = set.optimizer;
            obj.maxiter   = set.maxiter;
            obj.printing  = set.printing;
            obj.printMode = set.printMode;
            obj.createHistoyPrinter();
            
            obj.createMonitorDocker(cParams,set);
        end
        
        function itHas = hasFinished(obj,istep,nstep)
            itHas = obj.hasConverged || obj.niter >= obj.maxiter*(istep/nstep);
        end
        
    end
    
    methods (Access = private)
        
        function createHistoyPrinter(obj)
            settingsMetricsPrinter.fileName = obj.designVar.meshGiD.problemID;
            settingsMetricsPrinter.optimizer = obj;
            settingsMetricsPrinter.cost = obj.cost;
            settingsMetricsPrinter.constraint = obj.constraint;
            obj.historyPrinter = OptimizationMetricsPrinterFactory.create(obj.optimizer,obj.printing,settingsMetricsPrinter);
        end
        
        function createMonitorDocker(obj,cParams,set)
            mS.settingsParamsMonitor.showOptParams   = cParams.monitoring;
            mS.settingsParamsMonitor.convergenceVars = cParams.convergenceVars;
            mS.settingsParamsMonitor.settings        = set;
            mS.settingsDesignVarMonitor.plotting  = set.plotting;
            mS.settingsDesignVarMonitor.designVar = cParams.designVariable;
            mS.settingsDesignVarMonitor.settings  = set;
            
            obj.monitor = MonitoringDocker(mS);
        end
        
        function d = createPostProcessDataBase(obj,fileName)
            d.mesh    = obj.designVar.meshGiD;
            d.outName = fileName;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.getValue();
        end
        
    end
    
end
