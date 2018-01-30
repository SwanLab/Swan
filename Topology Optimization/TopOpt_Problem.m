classdef TopOpt_Problem < handle
    properties (GetAccess = public,SetAccess = public)
        cost
        constraint
        TOL
        x
        interpolation
        filter
        algorithm
        optimizer
        physicalProblem
        settings
        incropt
    end
    
    properties (Access = private)
        hole_value
        ini_design_value
    end
    
    methods (Static)
        function obj=create(settings)
            switch settings.ptype
                case 'Compliance_st_Volume'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_Compliance_st_Volume(settings);
                    obj.physicalProblem=Physical_Problem(obj.settings.filename);
                case 'Compliance_st_VolumePerimeter'
                    settings.nconstr=2;
                    obj=TopOpt_Problem_Compliance_st_VolumePerimeter(settings);
                    obj.physicalProblem=Physical_Problem(obj.settings.filename);
                case 'ComplianceLamPerimeter_st_Volume'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_ComplianceLamPerimeter_st_Volume(settings);
                    obj.physicalProblem=Physical_Problem(obj.settings.filename);
                case 'Chomog_alphabeta_st_Volume'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_Chomog(settings);
                    obj.physicalProblem=Physical_Problem_Micro(obj.settings.filename);
                case 'ChomogLamPerimeter_alphabeta_st_Volume'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_Chomog(settings);
                    obj.physicalProblem=Physical_Problem_Micro(obj.settings.filename);
                case 'Chomog_fraction_st_Volume'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_Chomog(settings);
                    obj.physicalProblem=Physical_Problem_Micro(obj.settings.filename);
                case 'ChomogLamPerimeter_fraction_st_Volume'
                    settings.nconstr=1;
                    obj=TopOpt_Problem_Chomog(settings);
                    obj.physicalProblem=Physical_Problem_Micro(obj.settings.filename);
                otherwise
                    error('Problem not added');
            end
        end
    end
    methods (Access = public)
        function obj=TopOpt_Problem(settings)
            obj.settings=settings;
            obj.TOL=obj.settings.TOL;
            obj.interpolation=Interpolation.create(obj.TOL,settings.material,settings.method);
            obj.ini_design_value=1;
            obj.hole_value=0;
            switch obj.settings.optimizer
                case 'SLERP'
                    obj.optimizer=Optimizer_AugLag(settings,Optimizer_SLERP(settings));
                    obj.ini_design_value=-1.015243959022692;
                    obj.hole_value=0.507621979511346;
                case 'PROJECTED GRADIENT'
                    obj.optimizer=Optimizer_AugLag(settings,Optimizer_PG(settings));
                case 'MMA'
                    obj.optimizer=Optimizer_MMA(settings);
                case 'IPOPT'
                    obj.optimizer=Optimizer_IPOPT(settings);
            end
            obj.filter=Filter.create(obj.settings.filter,obj.settings.optimizer);
        end
        
        function preProcess(obj)
            %initialize design variable
            obj.physicalProblem.preProcess;
            obj.filter.preProcess(obj.physicalProblem);
            obj.compute_initial_design;
            nsteps=obj.settings.nsteps;
            obj.incropt.alpha_vol = obj.generate_incr_sequence(1/nsteps,1,nsteps,'linear');
            obj.incropt.alpha_constr = obj.generate_incr_sequence(0,1,nsteps,'linear');
            obj.incropt.alpha_optimality= obj.generate_incr_sequence(0,1,nsteps,'linear');
        end
        function computeVariables(obj)
            for t = 1:obj.settings.nsteps
                incremental_step=t
                obj.update_target_parameters(t)
                obj.compute_physical_variables;
                obj.cost.computef(obj.x,obj.physicalProblem,obj.interpolation,obj.filter);
                obj.constraint.computef(obj.x, obj.physicalProblem, obj.interpolation,obj.filter);
                obj.optimizer.setPhysicalProblem(obj.physicalProblem);
                obj.x=obj.optimizer.solveProblem(obj.x,obj.cost,obj.constraint,obj.interpolation,obj.filter);
            end
        end
        function update_target_parameters(obj,t)
            target_parameters.Vfrac = (1-obj.incropt.alpha_vol(t))*obj.settings.Vfrac_initial+obj.incropt.alpha_vol(t)*obj.settings.Vfrac_final;
            % target_parameters.epsilon = update_parameter(target_parameters.epsilon_ini,target_parameters.epsilon_final,target_parameters.alpha_eps(t));
            %target_parameters.epsilon_isotropy = update_parameter(target_parameters.epsilon_isotropy_ini,target_parameters.epsilon_isotropy_final,target_parameters.alpha_isotropy2d(t));
            target_parameters.constr_tol = (1-obj.incropt.alpha_constr(t))*obj.settings.constr_initial+obj.incropt.alpha_constr(t)*obj.settings.constr_final;
            target_parameters.optimality_tol = (1-obj.incropt.alpha_optimality(t))*obj.settings.optimality_initial+obj.incropt.alpha_optimality(t)*obj.settings.optimality_final;
            obj.cost.target_parameters=target_parameters;
            %obj.cost.h_C_0=[];
            obj.constraint.target_parameters=target_parameters;
            obj.optimizer.target_parameters=target_parameters;
        end
        function postProcess(obj)
        end
    end
    methods (Static)
        function x = generate_incr_sequence (x1,x2,nsteps,type,factor)
            
            switch type
                case 'linear'
                    x = linspace(x1,x2,nsteps);
                    
                case 'epsilon_sequence'
                    frac = 2;
                    kmax = ceil(log10(x1/x2)/log10(frac));
                    x = epsilon0./frac.^(1:kmax);
                    
                case 'logarithmic'
                    x = logspace(log10(x1),log10(x2),nsteps);
                    
                case 'custom'
                    if nsteps < 2
                        x = x2;
                    else
                        isteps = 0:nsteps-1;
                        x = 1-(1-isteps/(nsteps-1)).^(factor);
                        x = (x2-x1)*x + x1;
                    end
                    
                case 'free'
                    x = zeros(1,nsteps);
                    x(end) = 1;
                    
                otherwise
                    error('Incremental sequence type not detected.')
            end
            
        end
    end
    
    methods (Access=private)
        function obj = compute_initial_design(obj)
            obj.x=obj.ini_design_value*ones(obj.physicalProblem.mesh.npnod,obj.physicalProblem.geometry.ngaus);
            switch obj.settings.initial_case
                case 'circle'
                    width = max(obj.physicalProblem.mesh.coord(:,1)) - min(obj.physicalProblem.mesh.coord(:,1));
                    height = max(obj.physicalProblem.mesh.coord(:,2)) - min(obj.physicalProblem.mesh.coord(:,2));
                    center_x = 0.5*(max(obj.physicalProblem.mesh.coord(:,1)) + min(obj.physicalProblem.mesh.coord(:,1)));
                    center_y = 0.5*(max(obj.physicalProblem.mesh.coord(:,2)) + min(obj.physicalProblem.mesh.coord(:,2)));
                    radius = 0.2*min([width,height]);
                    
                    initial_holes = (obj.physicalProblem.mesh.coord(:,1)-center_x).^2 + (obj.physicalProblem.mesh.coord(:,2)-center_y).^2 - radius^2 < 0;
                    obj.x(initial_holes) = obj.hole_value;
                    %fracc = 1;
                    
                case 'horizontal'
                    initial_holes = obj.physicalProblem.mesh.coord(:,2) > 0.6 | obj.physicalProblem.mesh.coord(:,2) < 0.4;
                    obj.x(initial_holes) = obj.hole_value;
                    %                   fracc = 1;
                case 'square'
                    width = max(obj.physicalProblem.mesh.coord(:,1)) - min(obj.physicalProblem.mesh.coord(:,1));
                    height = max(obj.physicalProblem.mesh.coord(:,2)) - min(obj.physicalProblem.mesh.coord(:,2));
                    center_x = 0.5*(max(obj.physicalProblem.mesh.coord(:,1)) + min(obj.physicalProblem.mesh.coord(:,1)));
                    center_y = 0.5*(max(obj.physicalProblem.mesh.coord(:,2)) + min(obj.physicalProblem.mesh.coord(:,2)));
                    
                    offset_x = 0.2*width;
                    offset_y = 0.2*height;
                    
                    xrange = obj.physicalProblem.mesh.coord(:,1) < (center_x+offset_x) & obj.physicalProblem.mesh.coord(:,1) > (center_x-offset_x);
                    yrange = obj.physicalProblem.mesh.coord(:,2) < (center_y+offset_y) & obj.physicalProblem.mesh.coord(:,2) > (center_y-offset_y);
                    initial_holes = and(xrange,yrange);
                    obj.x(initial_holes) = obj.hole_value;
                    %fracc = 1;
                    
                case 'feasible'
                    initial_holes = false(size(obj.physicalProblem.mesh.coord,1),1);
                    obj.x(initial_holes) = obj.hole_value;
                    %fracc = min(1,element.Vfrac);
                case 'rand'
                    initial_holes = rand(size(obj.physicalProblem.mesh.coord,1),1) > 0.1;
                    obj.x(initial_holes) = obj.hole_value;
                    %fracc = 1;
                case 'full'
                otherwise
                    error('Initialize design variable case not detected.');
            end
            rho_elem = obj.filter.getP0fromP1(obj.x);
            matprop = obj.interpolation.computeMatProp(rho_elem);
            obj.physicalProblem.setMatProps(matprop);
        end        
        function obj = compute_physical_variables(obj)
            switch obj.physicalProblem.mesh.scale
                case 'MICRO'
                    obj.physicalProblem.computeChomog;
                case 'MACRO'
                    obj.physicalProblem.computeVariables;
            end
        end
    end
end
