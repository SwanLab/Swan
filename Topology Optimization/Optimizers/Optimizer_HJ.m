classdef Optimizer_HJ < Optimizer_Unconstrained
    
    properties
        optimality_tol
        constr_tol
        HJiter
        HJiter0;
        HJiter_min = 1;
        niter = 1;
        e2
        % !! Move to ShFunc_Velocity (?) eventually !!
        filter
    end
    
    methods
        function obj = Optimizer_HJ(settings,epsilon)
            obj@Optimizer_Unconstrained(settings,epsilon);
            % !! Check wheter it affects the problem! !!
            %             obj.ini_design_value = -1.015243959022692;
            %             obj.hole_value = 0.507621979511346;
            obj.ini_design_value = -0.1;
            obj.hole_value = 0.1;
            
            obj.HJiter0 = settings.HJiter0;
            obj.HJiter = obj.HJiter0;
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
            x = obj.updatePhi(x_ini,obj.objfunc.gradient);
            cost.computef(x);
            constraint.computef(x);
            constraint = obj.setConstraint_case(constraint);
            obj.objfunc.computeFunction(cost,constraint)
            
            incr_norm_L2  = obj.norm_L2(x,x_ini);
            
            incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
            
            obj.stop_criteria = ~((incr_cost < 0 && incr_norm_L2 < obj.max_constr_change) || obj.kappa <= obj.kappa_min);
            
            obj.stop_vars(1,1) = incr_cost;     obj.stop_vars(1,2) = 0;
            obj.stop_vars(2,1) = incr_norm_L2;   obj.stop_vars(2,2) = obj.max_constr_change;
            obj.stop_vars(3,1) = obj.kappa;     obj.stop_vars(3,2) = obj.kappa_min;
            % !! CHANGE THIS !! Kappa is being taken as good even though
            % cost increases.
            
            if obj.stop_criteria
                obj.computeKappa;
            else
                obj.niter = obj.niter+1;
            end
        end
        
        function phi_vect = updatePhi(obj,design_variable,gradient)
            % !! PATCH !!
            load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));
            gradient = -obj.filter.regularize(design_variable,gradient);
            
            if length(dim) == 2
                load(fullfile(pwd,'Allaire_ShapeOpt','meshSize'));
                load(fullfile(pwd,'Allaire_ShapeOpt','RI'));
                
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
                [gradient_ALLAIRE, b_ALLAIRE] = regularize3(design_variable,gradient,dx,dy,dz);
                
                figure('NumberTitle', 'off', 'Name', 'ALLAIRE- b')
                subplot(2,2,1), surf(-b_ALLAIRE(:,:,2)), title('b_A - Root')
                subplot(2,2,2), surf(-b_ALLAIRE(:,:,end-1)), title('b_A - Tip')
                subplot(2,2,3), surf(permute(-b_ALLAIRE(ceil(size(b_ALLAIRE,1)/2),:,:),[2 3 1])), title('b_A - XY')
                subplot(2,2,4), surf(permute(-b_ALLAIRE(:,ceil(size(b_ALLAIRE,2)/2),:),[1 3 2])), title('b_A - XZ')
                
                for n = 1:length(design_variable)
                    phi(b1(n,1),b1(n,2),b1(n,3)) = design_variable(n);
                    V(b1(n,1),b1(n,2),b1(n,3)) = gradient(n);
                    V_ALLAIRE(b1(n,1),b1(n,2),b1(n,3)) = gradient_ALLAIRE(n);
                end
                
%                 figure('NumberTitle', 'off', 'Name', 'ALLAIRE- Regularized V')
%                 subplot(2,2,1), surf(V_ALLAIRE(:,:,2)), title('V_A - Root')
%                 subplot(2,2,2), surf(V_ALLAIRE(:,:,end-1)), title('V_A - Tip')
%                 subplot(2,2,3), surf(permute(V_ALLAIRE(ceil(size(V,1)/2),:,:),[2 3 1])), title('V_A - XY')
%                 subplot(2,2,4), surf(permute(V_ALLAIRE(:,ceil(size(V,2)/2),:),[1 3 2])), title('V_A - XZ')
%                 
%                 figure('NumberTitle', 'off', 'Name', 'FEM-MAT-OO - Regularized V')
%                 subplot(2,2,1), surf(V(:,:,2)), title('V - Root')
%                 subplot(2,2,2), surf(V(:,:,end-1)), title('V - Tip')
%                 subplot(2,2,3), surf(permute(V(ceil(size(V,1)/2),:,:),[2 3 1])), title('V - XY')
%                 subplot(2,2,4), surf(permute(V(:,ceil(size(V,2)/2),:),[1 3 2])), title('V - XZ')
%                 
%                 figure('NumberTitle', 'off', 'Name', 'FEM-MAT-OO/ALLAIRE')
%                 subplot(2,2,1), surf(V(:,:,2)./V_ALLAIRE(:,:,2)), title('V/V_A - Root')
%                 subplot(2,2,2), surf(V(:,:,end-1)./V_ALLAIRE(:,:,end-1)), title('V/V_A - Tip')
%                 subplot(2,2,3), surf(permute(V(ceil(size(V,1)/2),:,:)./V_ALLAIRE(ceil(size(V,1)/2),:,:),[2 3 1])), title('V/V_A - XY')
%                 subplot(2,2,4), surf(permute(V(:,ceil(size(V,2)/2),:)./V_ALLAIRE(:,ceil(size(V,2)/2),:),[1 3 2])), title('V/V_A - XZ')
                
%                 close; close; close;
%                 close;
                
                dt = 0.5*obj.e2*obj.kappa*min(dx,dy)/max(abs(gradient(:))) ;
                phi = solvelvlset3(phi,V,dt,obj.HJiter,0,30,5,dx,dy,dz);
                
                phi_vect(A1(:,:,:)) = phi(:,:,:);
                phi_vect = phi_vect';
            end
            
            % !! CHECK !!
            obj.opt_cond = obj.kappa;
        end
        
        function initKappa(obj,~,~,~)
            obj.kappa = 1;
            obj.HJiter = obj.HJiter0;
        end
        
        function computeKappa(obj)
            if obj.HJiter > obj.HJiter_min
                obj.HJiter = round(obj.HJiter/obj.kfrac);
            else
                obj.kappa = obj.kappa/obj.kfrac;
            end
        end
    end
end
