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
        
        function obj = TopOpt_Problem(settings)
            obj.createDesignVariable(settings);
            settings.pdim = obj.designVariable.mesh.pdim;
            obj.createHomogenizedVarComputer(settings)
            
            obj.createIncrementalScheme(settings);
            
            obj.cost = Cost(settings,obj.designVariable,obj.homogenizedVarComputer);
            obj.constraint = Constraint(settings,obj.designVariable,obj.homogenizedVarComputer);
            
            obj.createOptimizer(settings);
            
            obj.createVideoManager(settings);
        end
        
        
        function createOptimizer(obj,settings)
            obj.createOptimizerSettings(settings);
            obj.optimizer = Optimizer.create(obj.optimizerSettings);
        end
        
        function createOptimizerSettings(obj,settings)
            scS.filename        = settings.filename;
            scS.epsilon         = obj.incrementalScheme.targetParams.epsilon;
            scS.nVariables      = obj.designVariable.nVariables;
            
            lsS.scalarProductSettings = scS;
            lsS.line_search     = settings.line_search;
            lsS.optimizer       = settings.optimizer;
            lsS.HJiter0         = settings.HJiter0;
            lsS.filename        = settings.filename;
            lsS.kappaMultiplier = settings.kappaMultiplier;
            lsS.epsilon         = obj.incrementalScheme.targetParams.epsilon;
            
            
            uncOptimizerSettings = SettingsOptimizerUnconstrained();
            
            uncOptimizerSettings.lineSearchSettings    = lsS;
            uncOptimizerSettings.scalarProductSettings = scS;
            
            uncOptimizerSettings.e2                  = settings.e2;
            uncOptimizerSettings.filter              = settings.filter;
            uncOptimizerSettings.printChangingFilter = settings.printChangingFilter;
            uncOptimizerSettings.filename            = settings.filename;
            uncOptimizerSettings.ptype               = settings.ptype;
            uncOptimizerSettings.lb                  = settings.lb;
            uncOptimizerSettings.ub                  = settings.ub;
            
            uncOptimizerSettings.target_parameters     = obj.incrementalScheme.targetParams;
            uncOptimizerSettings.designVariable     = obj.designVariable;

            
            optSet.uncOptimizerSettings = uncOptimizerSettings;
            
            optSet.nconstr              = settings.nconstr;
            optSet.constraint_case      = settings.constraint_case;
            optSet.optimizer            = settings.optimizer;
            optSet.maxiter              = settings.maxiter;
            
            optSet.printing             = settings.printing;
            optSet.printMode            = settings.printMode;
            
            optSet.monitoring           = settings.monitoring;
            optSet.monitoring_interval  = settings.monitoring_interval;
            optSet.plotting             = settings.plotting;
            optSet.pdim                 = settings.pdim;
            optSet.showBC               = settings.showBC;
            
            optSet.settings   = settings;
            
            
            
            optSet.designVar            = obj.designVariable;
            optSet.target_parameters    = obj.incrementalScheme.targetParams;
            
            optSet.cost       = obj.cost;
            optSet.constraint = obj.constraint;
            
            obj.optimizerSettings = optSet;
            
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
            designVarSettings = SettingsDesignVariable();
            designVarSettings.mesh = mesh;
            designVarSettings.type = settings.designVariable;
            designVarSettings.levelSetCreatorSettings       = settings.levelSetDataBase;
            designVarSettings.levelSetCreatorSettings.ndim  = mesh.ndim;
            designVarSettings.levelSetCreatorSettings.coord = mesh.coord;
            designVarSettings.levelSetCreatorSettings.type = settings.initial_case;
            switch designVarSettings.levelSetCreatorSettings.type
                case 'holes'
                    designVarSettings.levelSetCreatorSettings.dirichlet = mesh.dirichlet;
                    designVarSettings.levelSetCreatorSettings.pointload = mesh.pointload;
            end
            obj.designVariable = DesignVariable.create(designVarSettings);
        end
        
        function createHomogenizedVarComputer(obj,settings)
            settings.nelem = size(obj.designVariable.mesh.connec,1);
            s.type                   = settings.homegenizedVariablesComputer;
            s.interpolation          = settings.materialInterpolation;
            s.dim                    = settings.pdim;
            s.typeOfMaterial         = settings.material;
            s.constitutiveProperties = settings.TOL;
            s.vademecumFileName      = settings.vademecumFileName;
            s.nelem                  = settings.nelem;
            obj.homogenizedVarComputer = HomogenizedVarComputer.create(s);
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
            
            settingsIncrementalScheme.mesh = obj.designVariable.mesh;
            
            obj.incrementalScheme = IncrementalScheme(settingsIncrementalScheme);
        end
        
        function solveCurrentProblem(obj)
            iStep = obj.incrementalScheme.iStep;
            nSteps = obj.incrementalScheme.nSteps;
            obj.optimizer.solveProblem(iStep,nSteps);
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