classdef TopOpt_Problem < handle
    properties (GetAccess = public,SetAccess = public)
        cost
        constraint
        x
        filter %% !! Remove when scalar product not done in Optimizer
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
            %% !! This should be done in settings class !!
            settings.nconstr = length(settings.constraint);
            obj.cost = Cost(settings,settings.weights);
            obj.constraint = Constraint(settings);
            
            % This PhysProb is only gonna be used by filters & incremental -> no need of specifying MICRO or MACRO
            % Consider turning it into a more generic class like FEM
            obj.topOpt_params = Physical_Problem(settings.filename);
            obj.settings = settings;
            switch obj.settings.optimizer
                case 'SLERP'
                    obj.optimizer = Optimizer_AugLag(settings,Optimizer_SLERP(settings));
                case 'PROJECTED GRADIENT'
                    obj.optimizer = Optimizer_AugLag(settings,Optimizer_PG(settings));
                case 'MMA'
                    obj.optimizer = Optimizer_MMA(settings);
                case 'IPOPT'
                    obj.optimizer = Optimizer_IPOPT(settings);
            end
            obj.filter = Filter.create(obj.settings.filter,obj.settings.optimizer);
            obj.optimizer.mesh=obj.topOpt_params.mesh; %% !!JUST TO MAKE PLOTTING WORK
        end
        
        function preProcess(obj)
            % Initialize design variable
            obj.topOpt_params.preProcess;
            obj.filters_preProcess;
            
            obj.incremental_scheme = Incremental_Scheme(obj.settings,obj.topOpt_params);
            obj.compute_initial_design;
            obj.topOpt_params = []; % !! To check that only it is used once (Debugging) !!
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
                files_name = obj.settings.filename;
                files_folder = fullfile(pwd,'Output');
                iterations = 0:obj.optimizer.niter;
                video_name = strcat('./Videos/Video_',obj.settings.ptype,'_',obj.settings.optimizer,'_',obj.settings.method,'_',int2str(obj.settings.nsteps) ...
                    ,'_0dot',int2str(10*obj.settings.Vfrac_final),'_',int2str(obj.optimizer.niter),'.gif');
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
            nukn = 1;
            dof_filter = DOF(obj.topOpt_params.problemID,obj.topOpt_params.geometry.nnode,obj.topOpt_params.mesh.connec,nukn,obj.topOpt_params.mesh.npnod,obj.topOpt_params.mesh.scale);
            switch obj.topOpt_params.mesh.scale
                case 'MACRO'
                    dof_filter.dirichlet = [];
                    dof_filter.dirichlet_values = [];
                    dof_filter.neumann = [];
                    dof_filter.neumann_values  = [];
                    dof_filter.constrained = dof_filter.compute_constrained_dof(obj.topOpt_params.mesh.scale);
                    dof_filter.free = dof_filter.compute_free_dof();
                case 'MICRO'
                    dof_filter.dirichlet = [];
                    dof_filter.dirichlet_values = [];
                    dof_filter.neumann = [];
                    dof_filter.neumann_values  = [];
                    dof_filter.constrained = dof_filter.compute_constrained_dof(obj.topOpt_params.mesh.scale);
                    dof_filter.free = dof_filter.compute_free_dof();
            end
            obj.topOpt_params.setDof(dof_filter)
            
            filter_params = obj.getFilterParams(obj.topOpt_params);
            obj.filter.preProcess(filter_params); % !! REMOVE WHEN SCALAR PRODUCT NOT IN OPTIMIZER !!
            obj.cost.preProcess(filter_params);
            obj.constraint.preProcess(filter_params);
        end
        
        function checkDerivative(obj)
            obj.preProcess;
            Msmooth = obj.filter.Msmooth;
            x0 = obj.x;
            % Initialize function
            epsi = 1e-6;
            %initial
            compliance0 = ShFunc_Compliance(obj.settings);
            % compliance0.h_C_0 = 1;
            volume0 = ShFunc_Volume(obj.settings);
            perimeter0 = ShFunc_Perimeter(obj.settings);
            %new
            compliance = ShFunc_Compliance(obj.settings);
            % compliance.h_C_0 = 1;
            volume = ShFunc_Volume(obj.settings);
            perimeter = ShFunc_Perimeter(obj.settings);
            
            obj.incremental_scheme.update_target_parameters(1,compliance0,volume0,perimeter0);
            obj.incremental_scheme.update_target_parameters(1,compliance,volume,perimeter);
            %evaluate initial
            compliance0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            volume0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            perimeter0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            
            compliance0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            volume0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            perimeter0.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            
            compliance.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            volume.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            perimeter.computef(x0,obj.topOpt_params,obj.interpolation,obj.filter);
            
            nnod = length(compliance0.gradient);
            g = zeros(nnod,1);
            gp = zeros(nnod,1);
            gv = zeros(nnod,1);
            for inode = 1:nnod
                if mod(inode,100) == 0
                    disp(strcat('Node: ',int2str(inode)));
                end
                xnew = x0;
                xnew(inode) = xnew(inode)-epsi;
                compliance.computef(xnew,obj.topOpt_params,obj.interpolation,obj.filter);
                volume.computef(xnew,obj.topOpt_params,obj.interpolation,obj.filter);
                perimeter.computef(xnew,obj.topOpt_params,obj.interpolation,obj.filter);
                g(inode) = (compliance0.value-compliance.value)/epsi;
                gv(inode) = (volume0.value-volume.value)/epsi;
                gp(inode) = (perimeter0.value-perimeter.value)/epsi;
            end
            fprintf('Relative error Volume: %g\n',obj.error_norm_field(gv,volume0.gradient,Msmooth));
            fprintf('Relative error Perimeter: %g\n',obj.error_norm_field(gp,perimeter0.gradient,Msmooth));
            fprintf('Relative error Compliance: %g\n',obj.error_norm_field(g,compliance0.gradient,Msmooth));
        end
        function enorm = error_norm_field(obj,gp,gp0,Msmooth)
            
            enodal = gp0 - gp;
            enorm = (enodal'*Msmooth*enodal)/(gp'*Msmooth*gp);
            
        end
    end
    
    methods (Access = private)
        function obj = compute_initial_design(obj)
            %% !! INCLUDE THIS INSIDE CLASS PHYSICAL_PROBLEM OR PARENT/CHILD!!
            switch obj.settings.optimizer
                case 'SLERP'
                    obj.ini_design_value = -2;
                    obj.hole_value = 0.507621979511346;
                otherwise
                    obj.ini_design_value = 1;
                    obj.hole_value = 0;
                    
            end
            obj.x = obj.ini_design_value*ones(obj.topOpt_params.mesh.npnod,1);
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
            obj.optimizer.Msmooth = obj.filter.Msmooth;
            obj.optimizer.Ksmooth = obj.filter.Ksmooth;
            obj.optimizer.epsilon_scalar_product_P1 = obj.incremental_scheme.epsilon;            
            %% !! COULD BE CLEANER, NOT IN IF --> O.O.P.  !!
            if strcmp(obj.settings.optimizer,'SLERP')
                sqrt_norma = obj.optimizer.scalar_product(obj.x,obj.x);
                obj.x = obj.x/sqrt(sqrt_norma);
            end
        end
    end
    methods (Access = private, Static)
        function filter_params = getFilterParams(topOpt_params)
            for igauss = 1:topOpt_params.geometry.ngaus
                filter_params.M0{igauss} = sparse(1:topOpt_params.mesh.nelem,1:topOpt_params.mesh.nelem,topOpt_params.geometry.dvolu(:,igauss));
            end
            filter_params.dof = topOpt_params.dof;
            filter_params.element = topOpt_params.element;
            filter_params.dvolu = sparse(1:topOpt_params.mesh.nelem,1:topOpt_params.mesh.nelem,sum(topOpt_params.geometry.dvolu,2));
            [filter_params.Ksmooth, filter_params.Msmooth] = topOpt_params.computeKM(2);
            filter_params.coordinates = topOpt_params.mesh.coord;
            filter_params.connectivities = topOpt_params.mesh.connec;
            filter_params.nelem = topOpt_params.mesh.nelem;
            filter_params.nnode = topOpt_params.geometry.nnode;
            filter_params.npnod = topOpt_params.mesh.npnod;
            filter_params.ngaus = topOpt_params.geometry.ngaus;
            filter_params.shape = topOpt_params.geometry.shape;
        end
    end
end

