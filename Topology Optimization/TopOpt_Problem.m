classdef TopOpt_Problem < handle
    properties (GetAccess = public,SetAccess = public)
        cost
        constraint
        x
        algorithm
        optimizer
        topOpt_params %% !! Only used for Ksmooth, Msmooth and geom/mesh/dofs (for filters & incremental)
        settings
        incremental_scheme
    end
    
    properties (Access = private)
        hole_value
        ini_design_value
    end
    
    methods (Access = public)
        function obj = TopOpt_Problem(settings)
            obj.cost = Cost(settings,settings.weights); % Change to just enter settings
            obj.constraint = Constraint(settings);
            
            % Consider turning it into a more generic class like FEM !!
                    obj.topOpt_params = DiffReact_Problem(settings.filename);
            obj.settings = settings;

            obj.incremental_scheme = Incremental_Scheme(obj.settings,obj.topOpt_params.mesh);
            switch obj.settings.optimizer
                case 'SLERP'
                    obj.optimizer = Optimizer_AugLag(settings,obj.topOpt_params.mesh,Optimizer_SLERP(settings,obj.incremental_scheme.epsilon));
                case 'PROJECTED GRADIENT'
                    obj.optimizer = Optimizer_AugLag(settings,obj.topOpt_params.mesh,Optimizer_PG(settings,obj.incremental_scheme.epsilon));
                case 'MMA'
                    obj.optimizer = Optimizer_MMA(settings,obj.topOpt_params.mesh);
                case 'IPOPT'
                    obj.optimizer = Optimizer_IPOPT(settings,obj.topOpt_params.mesh);
            end
        end
        
        function preProcess(obj)
            % Initialize design variable
            obj.topOpt_params.preProcess;
            obj.filters_preProcess;
            obj.compute_initial_design;
            obj.topOpt_params = []; %% !! To check that only it is used once (DEBUGGING) !!
        end
        
        function computeVariables(obj)
            for istep = 1:obj.settings.nsteps
                disp(strcat('Incremental step: ',int2str(istep)))
                obj.incremental_scheme.update_target_parameters(istep,obj.cost,obj.constraint,obj.optimizer);
                obj.cost.computef(obj.x);
                obj.constraint.computef(obj.x);
                obj.x = obj.optimizer.solveProblem(obj.x,obj.cost,obj.constraint);
            end
        end
        
        function postProcess(obj)
            % Video creation
            if obj.settings.printing
                gidPath = 'C:\Program Files\GiD\GiD 13.0.2';% 'C:\Program Files\GiD\GiD 13.0.3';
                files_name = obj.settings.case_file;
                files_folder = fullfile(pwd,'Output',obj.settings.case_file);
                iterations = 0:obj.optimizer.niter;
                video_name = strcat('./Videos/Video_',obj.settings.case_file,'_',int2str(obj.optimizer.niter),'.gif');
                My_VideoMaker = VideoMaker_TopOpt.Create(obj.settings.optimizer);
                My_VideoMaker.Set_up_make_video(gidPath,files_name,files_folder,iterations)
                %
                output_video_name_design_variable = fullfile(pwd,video_name);
                My_VideoMaker.Make_video_design_variable(output_video_name_design_variable)
                
                % %
                % output_video_name_design_variable_reg = fullfile(pwd,'DesignVariable_Reg_Video');
                % My_VideoMaker.Make_video_design_variable_reg(output_video_name_design_variable_reg)
                %
                % output_video_name_design_variable_reg = fullfile(pwd,'DesignVariable_Reg_Video');
                % My_VideoMaker.Make_video_design_variable_reg(output_video_name_design_variable_reg)
                %
                % output_video_name_stress = fullfile(pwd,'Stress_Video');
                % My_VideoMaker.Make_video_stress(output_video_name_stress)
            end
            
        end
        
        function obj = filters_preProcess(obj)
            obj.cost.preProcess;
            obj.constraint.preProcess;
        end       
    end
    
    methods (Access = private)
        function obj = compute_initial_design(obj)
            %% !! INCLUDE THIS INSIDE CLASS PHYSICAL_PROBLEM OR PARENT/CHILD!!
            switch obj.settings.optimizer
                case 'SLERP'
                    obj.ini_design_value = -1.015243959022692;
                    obj.hole_value = 0.507621979511346;
                otherwise
                    obj.ini_design_value = 1;
                    obj.hole_value = 0;
                    
            end
            obj.x = obj.ini_design_value*ones(obj.topOpt_params.geometry.interpolation.npnod,1);
            switch obj.settings.initial_case
                case 'circle'
                    width = max(obj.topOpt_params.mesh.coord(:,1)) - min(obj.topOpt_params.mesh.coord(:,1));
                    height = max(obj.topOpt_params.mesh.coord(:,2)) - min(obj.topOpt_params.mesh.coord(:,2));
                    center_x = 0.5*(max(obj.topOpt_params.mesh.coord(:,1)) + min(obj.topOpt_params.mesh.coord(:,1)));
                    center_y = 0.5*(max(obj.topOpt_params.mesh.coord(:,2)) + min(obj.topOpt_params.mesh.coord(:,2)));
                    radius = 0.2*min([width,height]);
                    
                    initial_holes = (obj.topOpt_params.mesh.coord(:,1)-center_x).^2 + (obj.topOpt_params.mesh.coord(:,2)-center_y).^2 - radius^2 < 0;
                    obj.x(initial_holes) = obj.hole_value;
                    %fracc = 1;
                    
                case 'horizontal'
                    initial_holes = obj.topOpt_params.mesh.coord(:,2) > 0.6 | obj.topOpt_params.mesh.coord(:,2) < 0.4;
                    obj.x(initial_holes) = obj.hole_value;
                    %                   fracc = 1;
                case 'square'
                    width = max(obj.topOpt_params.mesh.coord(:,1)) - min(obj.topOpt_params.mesh.coord(:,1));
                    height = max(obj.topOpt_params.mesh.coord(:,2)) - min(obj.topOpt_params.mesh.coord(:,2));
                    center_x = 0.5*(max(obj.topOpt_params.mesh.coord(:,1)) + min(obj.topOpt_params.mesh.coord(:,1)));
                    center_y = 0.5*(max(obj.topOpt_params.mesh.coord(:,2)) + min(obj.topOpt_params.mesh.coord(:,2)));
                    
                    offset_x = 0.2*width;
                    offset_y = 0.2*height;
                    
                    xrange = obj.topOpt_params.mesh.coord(:,1) < (center_x+offset_x) & obj.topOpt_params.mesh.coord(:,1) > (center_x-offset_x);
                    yrange = obj.topOpt_params.mesh.coord(:,2) < (center_y+offset_y) & obj.topOpt_params.mesh.coord(:,2) > (center_y-offset_y);
                    initial_holes = and(xrange,yrange);
                    obj.x(initial_holes) = obj.hole_value;
                    %fracc = 1;
                    
                case 'feasible'
                    initial_holes = false(size(obj.topOpt_params.mesh.coord,1),1);
                    obj.x(initial_holes) = obj.hole_value;
                    %fracc = min(1,element.Vfrac);
                case 'rand'
                    initial_holes = rand(size(obj.topOpt_params.mesh.coord,1),1) > 0.1;
                    obj.x(initial_holes) = obj.hole_value;
                    %fracc = 1;
                case 'full'
                otherwise
                    error('Initialize design variable case not detected.');
                    
            end
            %% PROVISIONAL 
            if strcmp(obj.settings.optimizer,'SLERP')
               sqrt_norma = obj.optimizer.optimizer_unconstr.scalar_product.computeSP(obj.x,obj.x);
               obj.x = obj.x/sqrt(sqrt_norma);
            end
        end
    end
end

