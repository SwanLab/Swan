classdef Optimizer_HJ < Optimizer_Unconstrained
    
    properties
        optimality_tol
        constr_tol
        HJiter
        HJiter_min = 1;
    end
    
    methods
        function obj = Optimizer_HJ(settings,epsilon)
            % !! PACHT !!
            HJiter0 = 30;
            obj@Optimizer_Unconstrained(settings,epsilon);
            %             obj.ini_design_value = -1.015243959022692;
            %             obj.hole_value = 0.507621979511346;
            % !! Currently NOT USED because init_design is loaded !!
            obj.ini_design_value = -0.1;
            obj.hole_value = 0.1;
            obj.HJiter = HJiter0;
            obj.kappa = 1;
            obj.kappa_min = 1e-5;
            obj.max_constr_change = +Inf;
            obj.kfrac = 2;
            obj.nconstr = settings.nconstr;
        end
        
        function optimality_tol = get.optimality_tol(obj)
            optimality_tol = (0.0175/1e-3)*obj.target_parameters.optimality_tol;
        end
        
        function constr_tol = get.constr_tol(obj)
            constr_tol(1:obj.nconstr) = obj.target_parameters.constr_tol;
        end
        
        function x = updateX(obj,x_ini,cost,constraint)
            % !! PATCH !!
            load(fullfile(pwd,'Allaire_ShapeOpt','meshSize'));
            load(fullfile(pwd,'Allaire_ShapeOpt','constraints_weights'));
            V = obj.objfunc.gradient/(dx*dy);
            V = obj.regularize(x_ini,V);
            
            dt = 0.5*obj.kappa*min(dx,dy)/max(abs(V(:))) ;
            
            x = obj.updatePhi(x_ini,V,dt);
            cost.computef(x);
            constraint.computef(x);
            constraint = obj.setConstraint_case(constraint);
            obj.objfunc.computeFunction(cost,constraint)
            
            incr_norm_L2  = obj.norm_L2(x,x_ini);
            incr_cost = (obj.objfunc.value - obj.objfunc.value_initial)/abs(obj.objfunc.value_initial);
            
            if obj.HJiter > obj.HJiter_min
                obj.HJiter = round(obj.HJiter/obj.kfrac);
            else
                obj.kappa = obj.kappa/obj.kfrac;
            end
            obj.stop_criteria = ~((incr_cost < 0 && incr_norm_L2 < obj.max_constr_change) || obj.kappa <= obj.kappa_min);
            
            obj.stop_vars(1,1) = incr_cost;     obj.stop_vars(1,2) = 0;
            obj.stop_vars(2,1) = incr_norm_L2;   obj.stop_vars(2,2) = obj.max_constr_change;
            obj.stop_vars(3,1) = obj.kappa;     obj.stop_vars(3,2) = obj.kappa_min;
        end
        
        function phi_vect = updatePhi(obj,design_variable,gradient,dt)
            % !! PATCH !!
            load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));
            load(fullfile(pwd,'Allaire_ShapeOpt','meshSize'));
            load(fullfile(pwd,'Allaire_ShapeOpt','constraints_weights'));
            load(fullfile(pwd,'Allaire_ShapeOpt','RI'));
            
            for n = 1:length(design_variable)
                phi(b1(n,1),b1(n,2)) = design_variable(n);
                V(b1(n,1),b1(n,2)) = gradient(n);
            end
            
            % !! Using Allaire's curvature instead of perimeter !!
            phi = solvelvlset(phi,V,dt,obj.HJiter,lagP,RIiter,RIfreq,dx,dy);
            phi_vect(A1(:,:)) = phi(:,:);
            phi_vect = phi_vect';
            
            % !! CHECK !!
            obj.opt_cond = obj.kappa;
        end
        
        function computeKappa(obj,~,~,~)
            obj.kappa = 1;
            % !! PATCH !!
            HJiter0 = 30;
            obj.HJiter = HJiter0;
        end
        
        function v = regularize(~,x,V_vect)
            % !! PATCH !!
            load(fullfile(pwd,'Allaire_ShapeOpt','conversion'));
            load(fullfile(pwd,'Allaire_ShapeOpt','meshSize'));
            
            for n = 1:length(x)
                phi(b1(n,1),b1(n,2)) = x(n);
            end
            
            for n = 1:length(V_vect)
                V(b1(n,1),b1(n,2)) = V_vect(n);
            end
            
            % Now we calculate the surface Dirac function.
            epsperim = min(dx,dy)/20 ;
            sx = phi./sqrt(phi.^2+epsperim^2) ;
            
            sxn = shift2n('n',sx,'ng') ;
            sxs = shift2n('s',sx,'ng') ;
            sxe = shift2n('e',sx,'ng') ;
            sxw = shift2n('w',sx,'ng') ;
            
            %We now calculate d(phi)/dx and d(phi)/dy:
            dsxx = (sxw-sxe)/(2*dx) ;
            dsxy = (sxn-sxs)/(2*dy) ;
            
            delta = 0.5*sqrt(dsxx.^2+dsxy.^2); % Surface Dirac function
            
            b = -V.*delta; % The right hand side is defined on the boundary so we use the delta function.
            
            b_vect(A1(:,:)) = b(:,:);
            b_vect = b_vect';
            
            filter = Filter_PDE_Density('Bridge_Quadrilateral','MACRO');
            filter.preProcess;
            filter.updateEpsilon(0.03);
            
            v = filter.getP1fromP1(b_vect);
        end
    end
end
