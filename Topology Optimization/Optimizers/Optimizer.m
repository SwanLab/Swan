classdef Optimizer < handle
    
    properties (Access = protected)
        designVariable
        dualVariable
        cost
        constraint
        outputFunction
        maxIter
        nIter = 0
        targetParameters
        dualUpdater
        primalUpdater
        constraintCase
        historyPrinter
        convergenceVars
        monitor
        postProcess
    end
    
    properties (GetAccess = public, SetAccess = protected, Abstract)
        type
    end
    
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            f   = OptimizerFactory();
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = protected)
        
        function initOptimizer(obj,cParams)
            obj.nIter             = 0;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.maxIter           = cParams.maxIter;
            obj.targetParameters  = cParams.targetParameters;
            obj.constraintCase    = cParams.constraintCase;
            obj.outputFunction    = cParams.outputFunction.monitoring;

            obj.createHistoryPrinter(cParams.historyPrinterSettings);
            obj.createConvergenceVariables(cParams);
            obj.createMonitorDocker(cParams.monitoringDockerSettings);
            obj.createPostProcess(cParams.postProcessSettings);
        end

        function createPrimalUpdater(obj,cParams)
            f                 = PrimalUpdaterFactory();
            obj.primalUpdater = f.create(cParams);
        end

        function createDualUpdater(obj,cParams)
            f               = DualUpdaterFactory();
            obj.dualUpdater = f.create(cParams);      
        end
        
        function isAcceptable = checkConstraint(obj)
            for i = 1:length(obj.constraint.value)
                switch obj.constraintCase{i}
                    case {'EQUALITY'}
                        fine = obj.checkEqualityConstraint(i);
                    case {'INEQUALITY'}
                        fine = obj.checkInequalityConstraint(i);
                end
                if fine
                    isAcceptable = true;
                else
                    isAcceptable = false;
                end
            end
        end

        function createHistoryPrinter(obj,cParams)
            cParams.optimizer  = obj;
            cParams.cost       = obj.cost;
            cParams.constraint = obj.constraint;
            cParams.optimizerName = cParams.optimizer.type;
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
            if cParams.shallPrint
                d = obj.createPostProcessDataBase(cParams);
                d.printMode = cParams.printMode;
                d.nDesignVariables = obj.designVariable.nVariables;
                obj.postProcess = Postprocess('TopOptProblem',d);
            end
        end

        function d = createPostProcessDataBase(obj,cParams)
            d.mesh    = obj.designVariable.mesh;
            d.outFileName = cParams.femFileName;
            d.ptype   = cParams.ptype;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.create();
            d.ndim       = obj.designVariable.mesh.ndim;
            d.pdim       = cParams.pdim;
            d.optimizer  = obj.type;
            d.cost       = obj.cost;
            d.constraint = obj.constraint;
            d.designVar  = obj.designVariable.type;
        end

    end

    methods (Access = private)

        function c = checkInequalityConstraint(obj,i)
            g = obj.constraint.value(i);
            c = g <= 0;
        end

        function c = checkEqualityConstraint(obj,i)
            g = obj.constraint.value(i);
            c = abs(g) < obj.targetParameters.constr_tol;
        end

    end

end