classdef TopOpt_Problem < handle
    properties (GetAccess = public,SetAccess = public)
        cost
        constraint
        x
        algorithm
        optimizer
        mesh
        settings
        incrementalScheme
        designVarInitializer
    end
    
    properties (Access = private)
        hole_value
        ini_design_value
    end
    
    methods (Access = public)
        function obj = TopOpt_Problem(settings)
            obj.mesh = Mesh_GiD(settings.filename);
            settings.pdim = obj.mesh.pdim;
            obj.settings = settings;
            obj.incrementalScheme = IncrementalScheme(settings,obj.mesh);
            obj.optimizer = OptimizerFactory().create(obj.settings.optimizer,settings,obj.mesh,obj.incrementalScheme.epsilon);
            obj.cost = Cost(settings,settings.weights);
            obj.constraint = Constraint(settings);
            obj.designVarInitializer = DesignVariableCreator(settings,obj.mesh);
        end
        
        function preProcess(obj)
            obj.cost.preProcess;
            obj.constraint.preProcess;
            obj.x = obj.designVarInitializer.getValue();
        end
        
        function computeVariables(obj)
            for istep = 1:obj.settings.nsteps
                obj.displayIncrementalIteration(istep)
                obj.incrementalScheme.update_target_parameters(istep,obj.cost,obj.constraint,obj.optimizer);
                obj.x = obj.optimizer.solveProblem(obj.x,obj.cost,obj.constraint,istep,obj.settings.nsteps);
            end
        end
        
        
        function postProcess(obj)
            % Video creation
            if obj.settings.printing
                gidPath = 'C:\Program Files\GiD\GiD 13.0.4';% 'C:\Program Files\GiD\GiD 13.0.3';
                files_name = obj.settings.case_file;
                files_folder = fullfile(pwd,'Output',obj.settings.case_file);
                iterations = 0:obj.optimizer.niter;
                video_name = strcat('./Videos/Video_',obj.settings.case_file,'_',int2str(obj.optimizer.niter),'.gif');
                My_VideoMaker = VideoMaker_TopOpt.Create(obj.settings.optimizer,obj.mesh.pdim,obj.settings.case_file);
                My_VideoMaker.Set_up_make_video(gidPath,files_name,files_folder,iterations)
                
                output_video_name_design_variable = fullfile(pwd,video_name);
                My_VideoMaker.Make_video_design_variable(output_video_name_design_variable)
            end
        end
    end
    
    methods (Access = private)
       
        function hasTo = hasToPrintIncrIter(obj)
            hasTo = obj.settings.printIncrementalIter;
            if isempty(obj.settings.printIncrementalIter)
                hasTo = true;
            end                            
        end
        
        function displayIncrementalIteration(obj,istep)
            if obj.hasToPrintIncrIter()
               disp(strcat('Incremental step: ',int2str(istep)))
            end                        
        end
        
    end
    
end