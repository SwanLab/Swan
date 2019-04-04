classdef TopOpt_Problem < handle
    
    properties (GetAccess = public, SetAccess = public)
        cost
        constraint
        designVariable
        x
        algorithm
        optimizer
        mesh
        settings
        incrementalScheme
    end
    
    properties (Access = private)
        hole_value
        ini_design_value
    end
    
    properties (Access = private)
        videoMaker
    end
        
    
    methods (Access = public)
        
        function obj = TopOpt_Problem(settings)
            obj.createDesignVariable(settings);
            
            settings.pdim = obj.mesh.pdim;
            obj.settings = settings;
            
            obj.createIncrementalScheme(settings);
            obj.optimizer = OptimizerFactory().create(settings.optimizer,settings,obj.designVariable,obj.incrementalScheme.targetParams.epsilon);
            obj.cost = Cost(settings,settings.weights);
            obj.constraint = Constraint(settings);
            
%             obj.createVideoMaker(settings);
        end
        
        function preProcess(obj)
            obj.cost.preProcess();
            obj.constraint.preProcess();
            obj.x = obj.designVariable.value;
        end
        
        function computeVariables(obj)
            obj.linkTargetParams();
            while obj.incrementalScheme.hasNext()
                obj.incrementalScheme.next();
                obj.solveCurrentProblem();
            end
        end
        
        
        function postProcess(obj)
            % Video creation
            if obj.settings.printing
                gidPath = 'C:\Program Files\GiD\GiD 13.0.4';% 'C:\Program Files\GiD\GiD 13.0.3';
                fileName = obj.settings.case_file;
                filePath = fullfile(pwd,'Output',obj.settings.case_file);
                iterations = 0:obj.optimizer.niter;
                videoName = strcat('./Videos/Video_',obj.settings.case_file,'_',int2str(obj.optimizer.niter),'.gif');
                
                videoMaker = VideoMakerTopOptFactory().create(obj.settings.case_file,obj.designVariable.type,obj.mesh.pdim);
                videoMaker.Set_up_make_video(gidPath,fileName,filePath,iterations)
                
                videoPath = fullfile(pwd,videoName);
                videoMaker.Make_video_design_variable(videoPath)
            end
        end
        
    end
    
    methods (Access = private)
        
        function createDesignVariable(obj,settings)
            obj.mesh = Mesh_GiD(settings.filename);
            designVarSettings.mesh = obj.mesh;
            designVarInitializer = DesignVariableCreator(settings,obj.mesh);
            designVarSettings.value = designVarInitializer.getValue();
            designVarSettings.optimizer = settings.optimizer;
            obj.designVariable = DesignVariableFactory().create(designVarSettings);
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
            
            settingsIncrementalScheme.mesh = obj.mesh;
            
            obj.incrementalScheme = IncrementalScheme(settingsIncrementalScheme);
        end
        
        function solveCurrentProblem(obj)
            istep = obj.incrementalScheme.iStep;
            obj.designVariable = obj.optimizer.solveProblem(obj.designVariable,obj.cost,obj.constraint,istep,obj.settings.nsteps);
            obj.x = obj.designVariable.value;
        end
        
        function linkTargetParams(obj)
            obj.cost.target_parameters = obj.incrementalScheme.targetParams;
            obj.constraint.target_parameters = obj.incrementalScheme.targetParams;
            obj.optimizer.target_parameters = obj.incrementalScheme.targetParams;
        end
%         
%         function createVideoMaker(obj,settings)
%             obj.videoMaker = VideoMakerFactory(settings.printing);                
%         end
        
    end
    
end