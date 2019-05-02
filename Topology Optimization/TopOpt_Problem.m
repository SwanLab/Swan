classdef TopOpt_Problem < handle
    
    properties (GetAccess = public, SetAccess = public)
        designVariable
        dualVariable
        cost
        constraint
        optimizer
        incrementalScheme
        optimizerSettings
    end
    
    properties (Access = private)
        videoManager
        homogenizedVarComputer
    end
    
    methods (Access = public)
        
        function obj = TopOpt_Problem(cParams)
            obj.createDesignVariable(cParams);
            obj.createDualVariable(cParams);            
            obj.createHomogenizedVarComputer(cParams)
            obj.createIncrementalScheme(cParams);
            
            obj.createCostAndConstraint(cParams);
            
            obj.createOptimizer(cParams);
            
            obj.createVideoManager(cParams);
        end
        
        
        
        function createOptimizer(obj,settings)
            obj.createOptimizerSettings(settings);
            obj.optimizer = Optimizer.create(obj.optimizerSettings);
        end
        
        function createOptimizerSettings(obj,cParams)
            s = cParams.optimizerSettings;
            s.uncOptimizerSettings.lineSearchSettings.scalarProductSettings.epsilon         = obj.incrementalScheme.targetParams.epsilon;
            s.uncOptimizerSettings.lineSearchSettings.scalarProductSettings.nVariables      = obj.designVariable.nVariables;
            s.uncOptimizerSettings.scalarProductSettings = s.uncOptimizerSettings.lineSearchSettings.scalarProductSettings;
            
            s.uncOptimizerSettings.lineSearchSettings.epsilon = obj.incrementalScheme.targetParams.epsilon;
            
            s.uncOptimizerSettings.targetParameters  = obj.incrementalScheme.targetParams;
            s.uncOptimizerSettings.designVariable     = obj.designVariable;
            
            s.designVar         = obj.designVariable;
            s.target_parameters = obj.incrementalScheme.targetParams;
            s.cost              = obj.cost;
            s.constraint        = obj.constraint;
            s.incrementalScheme = obj.incrementalScheme;
            s.dualVariable      = obj.dualVariable;
            
            obj.optimizerSettings = s;
            
        end
        
        function computeVariables(obj)
            while obj.incrementalScheme.hasNext()
                obj.incrementalScheme.next();
                obj.optimizer.solveProblem();
            end
        end
        
        function postProcess(obj)
            obj.videoManager.makeVideo(obj.optimizer.niter);
        end
        
    end
    
    methods (Access = private)
        
        function optSet = obtainOptimizersSettings(obj,settings)
            epsilon = obj.incrementalScheme.targetParams.epsilon;
            settings.optimizerSettings.uncOptimizerSettings.lineSearchSettings.epsilon = epsilon;
            settings.optimizerSettings.uncOptimizerSettings.scalarProductSettings.epsilon = epsilon;
            set = settings.clone();
            optSet = set.optimizerSettings;
        end
        
        function createDesignVariable(obj,cParams)
            s = cParams.designVarSettings;
            obj.designVariable = DesignVariable.create(s);
        end
        
        function createDualVariable(obj,cParams)
            cParamsD.nConstraints = numel(cParams.settings.constraint);
            obj.dualVariable = DualVariable(cParamsD);
        end        
        
        function createHomogenizedVarComputer(obj,cParams)
            s = cParams.homogenizedVarComputerSettings;
            obj.homogenizedVarComputer = HomogenizedVarComputer.create(s);
        end
        
        function createIncrementalScheme(obj,cParams)
            s = cParams.incrementalSchemeSettings;
            obj.incrementalScheme = IncrementalScheme(s);
        end
        
        function createCostAndConstraint(obj,set)
            cParams.settings = set.settings;
            cParams.designVariable = obj.designVariable;
            cParams.homogenizedVarComputer = obj.homogenizedVarComputer;
            cParams.targetParameters = obj.incrementalScheme.targetParams;
            obj.cost       = Cost(cParams);
            cParams.dualVariable = obj.dualVariable;
            obj.constraint = Constraint(cParams);
        end
        
        function createVideoManager(obj,cParams)
            s = cParams.videoManagerSettings;
            if s.shallPrint
                obj.videoManager = VideoManager(s);
            else
                obj.videoManager = VideoManager_Null(s);
            end
        end
        
    end
    
end
