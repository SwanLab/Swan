classdef TopOpt_Problem < handle

    properties (GetAccess = public, SetAccess = public)
        designVariable
        cost
        constraint
        optimizer
        algorithm
        incrementalScheme
        optimizerSettings
    end

    properties (Access = private)
        hole_value
        ini_design_value
        videoManager
        homogenizedVarComputer
    end

    methods (Access = public)

        function obj = TopOpt_Problem(cParams)
            obj.createDesignVariable(cParams);
            obj.createHomogenizedVarComputer(cParams)
            obj.createIncrementalScheme(cParams);

            obj.cost = Cost(cParams.settings,obj.designVariable,obj.homogenizedVarComputer);
            obj.constraint = Constraint(cParams.settings,obj.designVariable,obj.homogenizedVarComputer);

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

            s.uncOptimizerSettings.target_parameters  = obj.incrementalScheme.targetParams;
            s.uncOptimizerSettings.designVariable     = obj.designVariable;

            s.designVar         = obj.designVariable;
            s.target_parameters = obj.incrementalScheme.targetParams;
            s.cost              = obj.cost;
            s.constraint        = obj.constraint;

            obj.optimizerSettings = s;

        end


        function computeVariables(obj)
            obj.linkTargetParams();
            while obj.incrementalScheme.hasNext()
                obj.incrementalScheme.next();
                obj.solveCurrentProblem();
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

        function createHomogenizedVarComputer(obj,cParams)
            s = cParams.homogenizedVarComputerSettings;
            obj.homogenizedVarComputer = HomogenizedVarComputer.create(s);
        end

        function createIncrementalScheme(obj,cParams)
            s = cParams.incrementalSchemeSettings;
            obj.incrementalScheme = IncrementalScheme(s);
        end

        function solveCurrentProblem(obj)
            iStep  = obj.incrementalScheme.iStep;
            nSteps = obj.incrementalScheme.nSteps;
            obj.optimizer.solveProblem(iStep,nSteps);
        end

        function linkTargetParams(obj)
            obj.cost.target_parameters       = obj.incrementalScheme.targetParams;
            obj.constraint.target_parameters = obj.incrementalScheme.targetParams;
            obj.optimizer.target_parameters  = obj.incrementalScheme.targetParams;
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
