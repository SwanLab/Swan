classdef TopologyMonitoring < handle
    
    properties (Access = public)
        optimizerName
        designVariable
        cost
        constraint
        dualVariable
        historyPrinter
        convergenceVars
        monitorDocker
        postProcess
        incrementalScheme
        monitor
    end

    properties (Access = private)
        nIter
        hasFinished
        lineSearch
        lineSearchTrials
    end
    
    methods (Access = public)
        
        function obj = TopologyMonitoring(cParams)
            obj.init(cParams)         
        end
        
        function compute(obj,cParams)
            switch obj.optimizerName
                case 'fmincon'
                    obj.plotFmincon(cParams);
                case 'NullSpace'
                    obj.plotNullSpace(cParams);
                case 'AugmentedLagrangian'
                    obj.plotAugmentedLagrangian(cParams);
                case 'Bisection'
                    obj.plotBisection(cParams);
                case 'IPOPT'
                    obj.plotIPOPT(cParams);
                case 'MMA'
                    obj.plotMMA(cParams);
                otherwise
                    error('Optimizer not implemented')
            end
        end

        function create(obj,cParams)
            obj.createHistoryPrinter(cParams.historyPrinterSettings);
            obj.createConvergenceVariables(cParams);
            obj.createMonitorDocker(cParams.monitoringDockerSettings);
            obj.createPostProcess(cParams.postProcessSettings);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.optimizerName     = cParams.type;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVar;
            obj.dualVariable      = cParams.dualVariable;
            obj.incrementalScheme = cParams.incrementalScheme;
        end
        
        function plotFmincon(obj,cParams)
            obj.dualVariable.value = zeros(obj.constraint.nSF,1);
            obj.nIter              = cParams.nIter;
            obj.hasFinished        = cParams.hasFinished;
            foOpt                  = cParams.firstorderopt;
            normXsquare            = obj.designVariable.computeL2normIncrement();
            incX                   = sqrt(normXsquare);     
            obj.designVariable.updateOld();
            switch cParams.algorithm
                case 'sqp'
                    stepL = cParams.stepsize;
                case 'interior-point'
                    stepL = cParams.trustregionradius;
                otherwise
            end
            obj.printOptimizerVariable();
            obj.convergenceVars.reset();
            obj.convergenceVars.append(incX);
            obj.convergenceVars.append(foOpt);
            obj.convergenceVars.append(stepL);
            obj.refreshMonitoring();
            obj.printHistory();
        end
        
        function plotNullSpace(obj,cParams)
            deltaCost              = obj.cost.value - cParams.oldCost;
            obj.nIter              = cParams.nIter;
            normXsquare            = obj.designVariable.computeL2normIncrement();
            obj.lineSearch         = cParams.tau;
            obj.lineSearchTrials   = cParams.lineSearchTrials;
            obj.hasFinished        = cParams.hasFinished;
            incX                   = sqrt(normXsquare);
            obj.designVariable.updateOld();
            obj.printOptimizerVariable();
            obj.convergenceVars.reset();
            obj.convergenceVars.append(deltaCost);
            obj.convergenceVars.append(incX);
            obj.convergenceVars.append(obj.lineSearch);
            obj.convergenceVars.append(obj.lineSearchTrials);
            obj.refreshMonitoring();
            obj.printHistory();
        end
        
        function plotAugmentedLagrangian(obj,x,cParams)
            
        end
        
        function plotBisection(obj,x,cParams)
            
        end
        
        function plotIPOPT(obj,x,ccParams)
            
        end
        
        function plotMMA(obj,x,ccParams)
            obj.convergenceVars.reset();
            obj.convergenceVars.append(s.KKTnorm);
            obj.convergenceVars.append(s.outitFrac);
        end

        function createHistoryPrinter(obj,cParams)
            obj.historyPrinter = OptimizationMetricsPrinterFactory.create(cParams);
        end
        
        function createConvergenceVariables(obj,cParams)
            s                   = cParams.optimizerNames;
            cVarD               = ConvergenceVarsDispatcher.dispatchNames(s);
            n                   = numel(cVarD);
            cVar                = ConvergenceVariables(n);
            obj.convergenceVars = cVar;
        end
        
        function createMonitorDocker(obj,s)
            s.designVariable  = obj.designVariable;
            s.dualVariable    = obj.dualVariable;
            s.cost            = obj.cost;
            s.constraint      = obj.constraint;
            s.convergenceVars = obj.convergenceVars;
            obj.monitorDocker = MonitoringDocker(s);
        end
        
        function createPostProcess(obj,cParams)
            if cParams.shallPrint
                d                  = obj.createPostProcessDataBase(cParams);
                d.printMode        = cParams.printMode;
                d.nDesignVariables = obj.designVariable.nVariables;                
                obj.postProcess    = Postprocess('TopOptProblem',d);
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
               
        function printOptimizerVariable(obj)
            if ~isempty(obj.postProcess)
            d.fields     = obj.designVariable.getVariablesToPlot();
            d.cost       = obj.cost;
            d.constraint = obj.constraint;
            obj.postProcess.print(obj.nIter,d);
            end
        end
        
        function printHistory(obj)
            iStep = obj.incrementalScheme.iStep;
            obj.historyPrinter.print(obj.nIter,iStep);
        end
        
        function printHistoryFinalValues(obj)
            obj.historyPrinter.printFinal();
        end

        function refreshMonitoring(obj)
            iStep = obj.incrementalScheme.iStep;
            nStep = obj.incrementalScheme.nSteps;
            obj.monitorDocker.refresh(obj.nIter,obj.hasFinished,iStep,nStep);
        end
             
    end
    
end