classdef Optimizer_HJ < Optimizer_Unconstrained
    
    properties
        optimality_tol
        constr_tol
        HJiter
        HJiter0; % !! Could be set in settings !!
        HJiter_min = 1;
        % !! REMOVE?? !!
        allow
        niter_allow
        niter = 1;
        e2
        % !! Move to ShFunc_Velocity (?) eventually !!
        filter
        case_file %% !! PROVISIONAL: Just for 3D Shape Opt debugging !! Delete when done
    end
    
    methods
        function obj = Optimizer_HJ(settings,epsilon)
            obj@Optimizer_Unconstrained(settings,epsilon);
            % !! Check wheter it affects the problem! !!
            %             obj.ini_design_value = -1.015243959022692;
            %             obj.hole_value = 0.507621979511346;
            obj.ini_design_value = -0.1;
            obj.hole_value = 0.1;
            
            obj.case_file = settings.case_file;
            obj.HJiter0 = settings.HJiter0;
            obj.HJiter = obj.HJiter0;
            obj.allow = settings.allow;
            obj.niter_allow = settings.niter_allow;
            obj.e2 = settings.e2;
            obj.kappa = 1;
            obj.kappa_min = 1e-5;
            obj.max_constr_change = +Inf;
            obj.kfrac = 2;
            obj.nconstr = settings.nconstr;
            % !! Move to ShFunc_Velocity (?) eventually !!
            if strcmp(settings.filter,'P1')
                settings.filter = 'PDE';
                disp('Filter P1 changed to PDE for HJ velocity regularization');
            end
            obj.filter =  Filter.create(settings);
            obj.filter.preProcess;
            obj.filter.updateEpsilon(epsilon);
        end
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = (0.0175/1e-3)*obj.target_parameters.optimality_tol;
        end
        
        function constr_tol = get.constr_tol(obj)
            constr_tol(1:obj.nconstr) = obj.target_parameters.constr_tol;
        end
        
        function x = updateX(obj,x_ini,cost,constraint)
            if obj.niter > obj.niter_allow
                obj.allow = 0;
            end
            x = obj.updatePhi(x_ini,obj.objfunc.gradient);
            cost.computef(x);
            constraint.computef(x);
            constraint = obj.setConstraint_case(constraint);
            obj.objfunc.computeFunction(cost,constraint)
            
            incr_norm_L2  = obj.norm_L2(x,x_ini);
            
            %             incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
            incr_cost = (obj.objfunc.value - obj.objfunc.value_initial*(1+obj.allow))/abs(obj.objfunc.value_initial);
            
            obj.stop_criteria = ~((incr_cost < 0 && incr_norm_L2 < obj.max_constr_change) || obj.kappa <= obj.kappa_min);
            
            obj.stop_vars(1,1) = incr_cost;     obj.stop_vars(1,2) = 0;
            obj.stop_vars(2,1) = incr_norm_L2;   obj.stop_vars(2,2) = obj.max_constr_change;
            obj.stop_vars(3,1) = obj.kappa;     obj.stop_vars(3,2) = obj.kappa_min;
            
            if obj.stop_criteria
                if obj.HJiter > obj.HJiter_min
                    obj.HJiter = round(obj.HJiter/obj.kfrac);
                else
                    obj.kappa = obj.kappa/obj.kfrac;
                end
            else
                obj.niter = obj.niter+1;
            end
        end
        
        function phi_vect = updatePhi(obj,design_variable,gradient)
            % !! PATCH !!
            load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));
            
            if contains(lower(obj.case_file),'tri') || contains(lower(obj.case_file),'quad')
                load(fullfile(pwd,'Allaire_ShapeOpt','meshSize'));
                load(fullfile(pwd,'Allaire_ShapeOpt','RI'));
                gradient = -obj.filter.regularize(design_variable,gradient);
                
                for n = 1:length(design_variable)
                    phi(b1(n,1),b1(n,2)) = design_variable(n);
                    V(b1(n,1),b1(n,2)) = gradient(n);
                end
                
                dt = 0.5*obj.e2*obj.kappa*min(dx,dy)/max(abs(gradient(:))) ;
                
                % !! Using Allaire's curvature instead of perimeter !!
                phi = solvelvlset(phi,V,dt,obj.HJiter,0,RIiter,RIfreq,dx,dy);
                phi_vect(A1(:,:)) = phi(:,:);
                phi_vect = phi_vect';
                
                %             for n = 1:length(V)
                %                 V_mat(b1(n,1),b1(n,2)) = V(n);
                %             end
                %             figure, surf(V_mat);
            else
                dx = 1.25; dy = 1.25; dz = 1.25;
                gradient = regularize3(design_variable,gradient,dx,dy,dz);
                %                 gradient = -gradient;
                
                for n = 1:length(design_variable)
                    phi(b1(n,1),b1(n,2),b1(n,3)) = design_variable(n);
                    V(b1(n,1),b1(n,2),b1(n,3)) = gradient(n);
                end
                
                dt = 0.5*obj.e2*obj.kappa*min(dx,dy)/max(abs(gradient(:))) ;
                phi = solvelvlset3(phi,V,dt,obj.HJiter,0,30,5,dx,dy,dz);
                
                phi_vect(A1(:,:,:)) = phi(:,:,:);
                phi_vect = phi_vect';
            end
            
            % !! CHECK !!
            obj.opt_cond = obj.kappa;
        end
        
        function computeKappa(obj,~,~,~)
            obj.kappa = 1;
            obj.HJiter = obj.HJiter0;
        end
    end
end
