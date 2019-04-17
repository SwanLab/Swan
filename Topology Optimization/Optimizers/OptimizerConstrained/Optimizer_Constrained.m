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
        metricsPrinter
    end
    
    properties (Access = private)
        postProcess
        printing
        printMode
    end
    
    methods (Access = public)
        
        function obj = Optimizer_Constrained(cParams)
            set = cParams.settings.settings;
            designVar = cParams.designVariable;
            
            obj.constraintCase   = set.constraint_case;
            obj.hasConverged      = false;
            
            % !! REMOVE !!
            obj.cost = cParams.settings.cost;
            obj.constraint = cParams.settings.constraint;
            
            obj.init(set,designVar);
            
            mS.showOptParams = cParams.monitoring;
            mS.plotting  = set.plotting;
            mS.settings  = set;
            mS.designVar = cParams.designVariable;
            
            obj.monitor = MonitoringDocker(mS);
        end
        
        function designVar = solveProblem(obj,designVar,cost,constraint,istep,nstep)
            obj.createPostProcess(cost,constraint);
            x0 = designVar.value;
            cost.computeCostAndGradient(x0);
            constraint.computeCostAndGradient(x0);
            obj.print();
            
            obj.monitor.refresh(x0,obj.niter,cost,constraint,obj.stop_vars,obj.hasFinished(istep,nstep),istep,nstep);
            
            while ~obj.hasFinished(istep,nstep)
                obj.niter = obj.niter+1;
                x = obj.update(x0);
                designVar.update(x);
                obj.monitor.refresh(x,obj.niter,cost,constraint,obj.stop_vars,obj.hasFinished(istep,nstep),istep,nstep);
                obj.print();
                obj.exportMetrics(istep)
                x0 = x;
            end
            obj.print();
            
            obj.hasConverged = 0;
        end
        
    end
    
    methods (Access = protected)
        
        function print(obj)
            d.x = obj.designVar.value;
            d.cost = obj.cost;
            d.constraint = obj.constraint;
            obj.postProcess.print(obj.niter,d);
        end
        
        function exportMetrics(obj,iStep)
            obj.metricsPrinter.print(obj.niter,iStep);
        end
        
        function createPostProcess(obj,cost,constraint)
            fileName = obj.designVar.meshGiD.problemID;
            d = obj.createPostProcessDataBase(fileName);
            d.printMode = obj.printMode;
            d.optimizer = obj.optimizer;
            d.cost = cost;
            d.constraint = constraint;
            if obj.printing
                obj.postProcess = Postprocess('TopOptProblem',d);
            else
                obj.postProcess = Postprocess_Null('',d);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function itHas = hasFinished(obj,istep,nstep)
            itHas = obj.hasConverged || obj.niter >= obj.maxiter*(istep/nstep);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,settings,designVar)
            obj.designVar = designVar;
            obj.optimizer = settings.optimizer;
            obj.maxiter   = settings.maxiter;
            obj.printing  = settings.printing;
            obj.printMode = settings.printMode;
            obj.createMetricsPrinter();
        end
        
        function createMetricsPrinter(obj)
            settingsMetricsPrinter.fileName = obj.designVar.meshGiD.problemID;
            settingsMetricsPrinter.optimizer = obj;
            settingsMetricsPrinter.cost = obj.cost;
            settingsMetricsPrinter.constraint = obj.constraint;
            obj.metricsPrinter = OptimizationMetricsPrinterFactory.create(obj.optimizer,obj.printing,settingsMetricsPrinter);
        end
        
        function d = createPostProcessDataBase(obj,fileName)
            d.mesh    = obj.designVar.meshGiD;
            d.outName = fileName;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.getValue();
        end
        
    end
    
end
