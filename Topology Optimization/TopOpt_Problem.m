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
    end
    
    properties (Access = private)
        videoManager
    end
        
    
    methods (Access = public)
        
        function obj = TopOpt_Problem(settings)
            obj.createDesignVariable(settings);
            
            settings.pdim = obj.designVariable.meshGiD.pdim;
            
            obj.createIncrementalScheme(settings);

            
            obj.cost = Cost(settings,settings.weights,obj.designVariable);
            obj.constraint = Constraint(settings,obj.designVariable);            

            obj.createOptimizerSettings(settings);             
            obj.optimizer = OptimizerFactory.create(obj.optimizerSettings);

            obj.createVideoManager(settings);
        end
        
        function createOptimizerSettings(obj,settings)
            lsS.line_search     = settings.line_search;
            lsS.optimizer       = settings.optimizer;
            lsS.HJiter0         = settings.HJiter0;
            lsS.filename        = settings.filename;
            lsS.kappaMultiplier = settings.kappaMultiplier;
            lsS.epsilon         = obj.incrementalScheme.targetParams.epsilon;
            
            
            scS.filename        = settings.filename;
            scS.epsilon         = obj.incrementalScheme.targetParams.epsilon;
            
            uncOptimizerSettings = SettingsOptimizerUnconstrained();
            
            uncOptimizerSettings.target_parameters     = settings.target_parameters;
            uncOptimizerSettings.lineSearchSettings    = lsS;
            uncOptimizerSettings.scalarProductSettings = scS;
            
            uncOptimizerSettings.e2                  = settings.e2;
            uncOptimizerSettings.filter              = settings.filter;
            uncOptimizerSettings.printChangingFilter = settings.printChangingFilter;
            uncOptimizerSettings.filename            = settings.filename;
            uncOptimizerSettings.ptype               = settings.ptype;
            
            optSet.uncOptimizerSettings = uncOptimizerSettings;
            optSet.monitoring           = settings.monitoring;
            optSet.nconstr              = settings.nconstr;
            optSet.target_parameters    = settings.target_parameters;
            optSet.constraint_case      = settings.constraint_case;   
            optSet.optimizer            = settings.optimizer;
            optSet.maxiter              = settings.maxiter;
            optSet.printing             = settings.printing;
            optSet.printMode            = settings.printMode;            
            
            optSet.plotting             = settings.plotting;
            optSet.pdim                 = settings.pdim;
            optSet.showBC               = settings.showBC;   
            
            optSet.settings = settings;
            optSet.cost     = obj.cost;
            optSet.constraint = obj.constraint;
            
            optSet.designVar = obj.designVariable;
            obj.optimizerSettings = optSet;
            
        end
        
        function preProcess(obj)
            obj.cost.preProcess();
            obj.constraint.preProcess();
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
        
        function createDesignVariable(obj,settings)
            mesh = Mesh_GiD(settings.filename);
            designVarSettings.mesh = mesh;
            designVarInitializer = DesignVariableCreator(settings,mesh);
            designVarSettings.value = designVarInitializer.getValue();
            %designVarSettings.optimizer = settings.optimizer;
            designVarSettings.type = settings.designVariable;
            obj.designVariable = DesignVariable.create(designVarSettings);
        end
        
        function createIncrementalScheme(obj,settings)
            settingsIncrementalScheme = SettingsIncrementalScheme();
            
            settingsIncrementalScheme.settingsTargetParams.VfracInitial = settings.Vfrac_initial;
            settingsIncrementalScheme.settingsTargetParams.VfracFinal = settings.Vfrac_final;
            settingsIncrementalScheme.settingsTargetParams.constrInitial = settings.constr_initial;
            settingsIncrementalScheme.settingsTargetParams.constrFinal = settings.constr_final;
            settingsIncrementalScheme.settingsTargetParams.optimalityInitial = settings.optimality_initial;
            settingsIncrementalScheme.settingsTargetParams.optimalityFinal = settings.optimality_final;
            settingsIncrementalScheme.settingsTargetParams.epsilonInitial = settings.epsilon_initial;
            settingsIncrementalScheme.settingsTargetParams.epsilonFinal = settings.epsilon_final;
            settingsIncrementalScheme.settingsTargetParams.epsilonIsotropyInitial = settings.epsilon_isotropy_initial;
            settingsIncrementalScheme.settingsTargetParams.epsilonIsotropyFinal = settings.epsilon_isotropy_final;
            
            settingsIncrementalScheme.nSteps = settings.nsteps;
            settingsIncrementalScheme.shallPrintIncremental = settings.printIncrementalIter;
            
            settingsIncrementalScheme.mesh = obj.designVariable.meshGiD;
            
            obj.incrementalScheme = IncrementalScheme(settingsIncrementalScheme);
        end
        
        function solveCurrentProblem(obj)
            iStep = obj.incrementalScheme.iStep;
            nSteps = obj.incrementalScheme.nSteps;
            obj.designVariable = obj.optimizer.solveProblem(obj.designVariable,obj.cost,obj.constraint,iStep,nSteps);
        end
        
        function linkTargetParams(obj)
            obj.cost.target_parameters = obj.incrementalScheme.targetParams;
            obj.constraint.target_parameters = obj.incrementalScheme.targetParams;
            obj.optimizer.target_parameters = obj.incrementalScheme.targetParams;
        end
       
        function createVideoManager(obj,settings)
            if settings.printing
                obj.videoManager = VideoManager(settings,obj.designVariable);
            else
                obj.videoManager = VideoManager_Null(settings,obj.designVariable);
            end
        end
        
    end
    
end