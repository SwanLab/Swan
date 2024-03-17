classdef TopologyMonitoring < handle
    
    properties (Access = public)
        type
        designVariable
        cost
        constraint
        dualVariable
        historyPrinter
        convergenceVars
        monitorDocker
        postProcess
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
            switch obj.type
                case 'fmincon'
                    obj.plotFmincon(cParams);
                case 'NullSpace'
                    obj.plotNullSpace(cParams);
                case 'AlternatingPrimalDual'
                    obj.plotAugmentedLagrangian(cParams);
                case 'DualNestedInPrimal'
                    obj.plotBisection(cParams);
                case 'IPOPT'
                    obj.plotIPOPT(cParams);
                case 'MMA'
                    obj.plotMMA(cParams);
                case 'IPM'
                    obj.plotIPM(cParams);
                otherwise
                    error('Optimizer not implemented')
            end
        end

        function create(obj,cParams)
            %obj.createHistoryPrinter(cParams.historyPrinterSettings);
            %obj.createConvergenceVariables(cParams);
            %obj.createMonitorDocker(cParams.monitoringDockerSettings);
            %obj.createPostProcess(cParams.postProcessSettings);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.type     = cParams.type;
            obj.cost              = cParams.cost;
            obj.constraint        = cParams.constraint;
            obj.designVariable    = cParams.designVariable;
            obj.dualVariable      = cParams.dualVariable;
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
            obj.convergenceVars.append(cParams.meritNew);
            obj.refreshMonitoring();
            obj.printHistory();
        end
        
        function plotAugmentedLagrangian(obj,cParams)
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
            obj.convergenceVars.append(cParams.meritNew);
            obj.refreshMonitoring();
            obj.printHistory();
        end
        
        function plotBisection(obj,cParams)
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
        
        function plotIPOPT(obj,cParams)
            obj.hasFinished = cParams.hasFinished;
            obj.nIter       = cParams.nIter;
            normXsquare     = obj.designVariable.computeL2normIncrement();            
            obj.printOptimizerVariable();
            obj.convergenceVars.reset();
            obj.convergenceVars.append(cParams.inf_pr);
            obj.convergenceVars.append(cParams.inf_du);            
            obj.convergenceVars.append(sqrt(normXsquare));
            obj.refreshMonitoring();
            obj.printHistory();          
        end
        
        function plotMMA(obj,cParams)
            obj.hasFinished = cParams.hasFinished;
            obj.nIter       = cParams.nIter;
            obj.printOptimizerVariable();
            obj.convergenceVars.reset();
            obj.convergenceVars.append(cParams.KKTnorm);
            obj.convergenceVars.append(cParams.outitFrac);
            obj.refreshMonitoring();
            obj.printHistory();
        end

        function plotIPM(obj,cParams)
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
            obj.convergenceVars.append(cParams.meritNew);
            obj.refreshMonitoring();
            obj.printHistory();
        end

        function createHistoryPrinter(obj,cParams)
            cParams.optimizer  = obj;
            cParams.cost       = obj.cost;
            cParams.constraint = obj.constraint;
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
            d.ptype   = cParams.ptype;
            ps = PostProcessDataBaseCreator(d);
            d  = ps.create();
            d.pdim    = cParams.pdim;     
            d.ndim    = obj.designVariable.mesh.ndim;
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
%             obj.postProcess.print(obj.nIter,d);
            end
        end
        
        function printHistory(obj)
            obj.historyPrinter.print(obj.nIter,1);
        end
        
        function printHistoryFinalValues(obj)
            obj.historyPrinter.printFinal();
        end

        function refreshMonitoring(obj)
            obj.monitorDocker.refresh(obj.nIter,obj.hasFinished,1,1);
        end
             
    end
    
end